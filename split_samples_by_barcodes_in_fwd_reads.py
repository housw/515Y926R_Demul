#!/usr/bin/env python

# Copyright (C) 2019  Shengwei Hou
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


import os
import re
import sys
import gzip
import logging
import argparse
#import coloredlogs
import numpy as np
from collections import OrderedDict


# logging
logfmt = "[%(asctime)s] [%(levelname)s] [%(name)s] %(message)s"
#field_styles= {'asctime': {'color': 'green'}, 'levelname': {'bold': True, 'color': 'cyan'},
#               'name': {'color': 'blue'}, 'username': {'color': 'magenta'}, 'hostname': {'color': 'yellow'}}
datefmt = "%Y-%m-%d %H:%M:%S"
logging.basicConfig(stream=sys.stdout, format=logfmt, datefmt=datefmt)
logger = logging.getLogger('main')
#coloredlogs.install(fmt=logfmt, datefmt=datefmt, field_styles=field_styles)


def set_loglevel(loglevel):
    """
    :param loglevel: loglevel (str): minimum loglevel for emitting messages
    :return:
    """
    _level = {
        'critical': logging.CRITICAL,
        'error': logging.ERROR,
        'warning': logging.WARNING,
        'info': logging.INFO,
        'debug': logging.DEBUG}.get(loglevel, logging.DEBUG)

    logging.getLogger('main').setLevel(_level)



class Fastq(object):

    def __init__(self, header, seq, qual):
        self.header = header
        self.seq = seq.upper()
        self.qual = qual
        self.name = self.header.strip().split()[0]

    def print_pretty_fasta(self, width=100):
        """
        :param width: maximum 'width' per line for sequences
        :return:      return pretty string representation of Fasta record
        """
        ret = ">" + self.header + "\n"
        for i, char in enumerate(self.seq):
            if i > 0 and i % width == 0:
                ret += "\n"
            ret += char
        return ret + "\n"

    def __repr__(self):
        return 'Fastq(header=%s, seq=%s, qual=%s)' % (self.header, self.seq, self.qual)

    def __str__(self):
        return "@" + self.header + "\n" \
               + self.seq +"\n" \
               + "+\n" \
               + str(self.qual) + "\n"


def parse_fastq(fastq_file):
    """
    :param fastq_file: input fastq file
    :return:           yield Fastq record as a generator
    """

    # handle gzip compressed file
    if fastq_file.endswith(".gz"):
        ih = gzip.open(fastq_file, "rt")
    else:
        ih = open(fastq_file, "r")

    # parse fastq records
    header = ih.readline().strip()[1:]
    while header:
        header = header.strip()[1:]
        seq = ih.readline().strip()
        ih.readline()
        qual = ih.readline().strip()
        yield Fastq(header, seq, qual)
        header = ih.readline()

    ih.close()


# https://en.wikipedia.org/wiki/Hamming_distance#Algorithm_example
def hamming_distance(s1, s2):
    #Return the Hamming distance between equal-length sequences
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


def get_barcode2samples(barcode_file):
    """ read in two column barcode file, return barcode2sample dict
    """

    barcode2samples = {} # {barcode:sample}

    with open(barcode_file, "r") as ih:
        for oline in ih:
            if oline.startswith("#"):
                continue
            line = oline.strip().split()
            if len(line) < 2:
                print("Either sample or barcode is missing for {0}".format(oline))
                continue
            sample, barcode = line[0], line[1]
            assert barcode not in barcode2samples, "Barcode %s are not unique!"%barcode
            barcode2samples[barcode] = sample

    return barcode2samples                        


def get_barcode2handles(barcode2samples, prefix):
    """ get barcode to fwd_output and rev_output handles 
    """
    barcode2handles = {} #{barcode:[fwd_handle, rev_handle]}

    for barcode, sample in barcode2samples.items():
        fwd_handle = open(prefix+sample+".R1.fastq", "w")
        rev_handle = open(prefix+sample+".R2.fastq", "w") 
        barcode2handles[barcode] = [fwd_handle, rev_handle]
    
    # add unknown
    barcode2handles["unknown"] = [open(prefix+"unknown.R1.fastq", "w"), open(prefix+"unknown.R1.fastq", "w")]    

    return barcode2handles


def get_dist_1_close_barcode(barcode, barcodes):
    dist = []
    for bar in barcodes:
        #sys.stdout.write(barcode +"\t"+ ind +"\n")
        dist.append(hamming_distance(barcode, bar))
    dist = np.array(dist)
    # close barcode should be within 1 hamming distance, and should be the only best one
    if dist.min() == 1 and sum(dist == dist.min()) == 1:
        idx = np.where(dist==dist.min())[0][0]
        #print("idx is ", idx)
        close_barcode = barcodes[idx]
    else:
        close_barcode = "unknown"

    return close_barcode
    

def split_samples_by_barcode(fwd_fastq, rev_fastq, barcode2samples, prefix, start_position, barcode_length, primer_regex):
    """
    """

    barcode2handles = get_barcode2handles(barcode2samples, prefix)
    sample2counts = {}  # calcualte read count assigned to each sample
    barcode2counts = {} # calcualte read count assigned to detected barcodes 
    for sample in barcode2samples.values():
        sample2counts[sample] = 0
    sample2counts['unknown'] = 0

    for fwd_read, rev_read in zip(parse_fastq(fwd_fastq), parse_fastq(rev_fastq)):
        if fwd_read.name == rev_read.name:
            logger.debug("-----------------------------------")

            # parse barcode by primer regex matching
            match = re.search(primer_regex, fwd_read.seq)
            if match:
                primer_start = match.start()
                primer_seq = match.group()
                if primer_start+1 > barcode_length:
                    barcode_start = primer_start - barcode_length
                    barcode = fwd_read.seq[barcode_start:primer_start]
                else:
                    logger.warning("barcode sequence is truncated since primer "\
                                   "starts at 1-based index {0}".format(primer_start+1))
                    barcode = "unknown"
                logger.debug("barcode: {0}{1}".format("".join([" "]*barcode_start) , barcode))
                logger.debug("primer : {0}{1}".format("".join([" "]*primer_start), primer_seq))
                logger.debug("fwd seq: {0}".format(fwd_read.seq[:50]))
            else:
                barcode = fwd_read.seq[start_position-1:start_position-1+barcode_length]
                logger.debug("barcode: {0}{1}".format("".join([" "]*(start_position-1)), barcode))
                logger.debug("fwd seq: {0}".format(fwd_read.seq[:50]))
            

            if barcode in barcode2counts:
                barcode2counts[barcode] += 1
            else:
                barcode2counts[barcode] = 1
            if barcode in barcode2samples:
                barcode2handles[barcode][0].write(str(fwd_read))
                barcode2handles[barcode][1].write(str(rev_read))
                sample2counts[barcode2samples[barcode]] += 1 
            else:
                # find closest barcode
                if barcode != "unknown":
                    logger.debug("try to find the closest known barcode...")
                    close_barcode = get_dist_1_close_barcode(barcode, list(barcode2samples.keys()))
                    if close_barcode != "unknown":
                        logger.debug("found  : {0}{1}".format("".join([" "]*barcode_start) , close_barcode))
                        barcode2handles[close_barcode][0].write(str(fwd_read))
                        barcode2handles[close_barcode][1].write(str(rev_read))
                        sample2counts[barcode2samples[close_barcode]] += 1
                    else:
                        logger.debug("no close barcode found, set to unknown")
                        barcode = "unknown"
                if barcode == "unknown":
                    barcode2handles['unknown'][0].write(str(fwd_read))
                    barcode2handles['unknown'][1].write(str(rev_read))
                    sample2counts['unknown'] += 1
        else:
            logger.error("Forward read %s is different with reverse read %s"%(fwd_read.name, rev_read.name)) 
            barcode2handles['unknown'][0].write(str(fwd_read))
            barcode2handles['unknown'][1].write(str(rev_read))
            sample2counts['unknown'] += 1

    for handles in barcode2handles.values():
        for handle in handles:
            handle.close()
    with open(prefix+"sample2count.tsv", "w") as oh:
        oh.write("Sample\tCount\n")
        for sample in sorted(sample2counts.keys()):
            count = sample2counts[sample]
            oh.write(str(sample) +"\t"+ str(count) +"\n")
    with open(prefix+"barcode2count.tsv", "w") as oh:
        oh.write("Barcode\tCount\n")
        sorted_barcode2counts = OrderedDict(sorted(barcode2counts.items(), key=lambda x: x[1], reverse=True))
        for barcode, count in sorted_barcode2counts.items():
            oh.write(str(barcode) +"\t"+ str(count) +"\n")


def main():

    # main parser
    parser = argparse.ArgumentParser(description="split samples in fastq file based on barcodes found in forward reads at certain position") 
    parser.add_argument("barcode_file", help="input barcode file")
    parser.add_argument("input_fwd_fastq", help="input forward reads in fastq format")
    parser.add_argument("input_rev_fastq", help="input reverse reads in fastq format")
    parser.add_argument("-s", "--start_position", type=int, default=5, help="start position of barcodes in sequences, 1-based coordinate, default=5")
    parser.add_argument("-r", "--primer_regex", type=str, default="GTG[TC]CAGC[AC]GCCGCGGTAA", help="primer regex, default=GTG[TC]CAGC[AC]GCCGCGGTAA")
    parser.add_argument("-b", "--barcode_length", type=int, default=5, help="barcode length, default=5")
    parser.add_argument("-o", "--output_dir", default="./", help="output directory, default=./")
    parser.add_argument("-p", "--prefix", default="split_", help="output prefix, default=split_")
    parser.add_argument('-l', '--loglevel', default='info', choices=['critical', 'error', 'warning', 'info', 'debug'],
                        help='logging level')
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 1.0")


    if len(sys.argv) < 2:
        print("\nERROR: Not enough parameters were provided, please refer to the usage.\n", file=sys.stderr)
        print(parser.format_help(), file=sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    # logging
    set_loglevel(args.loglevel)

    # output dir handeling
    if not os.path.exists(args.output_dir):
        os.system("mkdir -p {}".format(args.output_dir))
    prefix = os.path.join(args.output_dir, args.prefix)

    # get barcode2samples
    logger.info("parsing barcode file ...")
    barcode2samples = get_barcode2samples(args.barcode_file)

    # split samples
    logger.info("\nSpliting {0} and {1} based on {2} ...".format(args.input_fwd_fastq, args.input_rev_fastq, args.barcode_file))
    split_samples_by_barcode(args.input_fwd_fastq, args.input_rev_fastq, barcode2samples, prefix, args.start_position, args.barcode_length, args.primer_regex)
    logger.info("Done!\n")


if __name__ == "__main__":
    main()
