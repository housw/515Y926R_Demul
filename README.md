# 515Y926R_Demul
A simple barcode splitter for 515F/926R amplified amplicon sequencing data. It assumes a pair of `fastq` files that are amplified using the Fuhrmam Lab [protocol](https://www.protocols.io/view/fuhrman-lab-515f-926r-16s-and-18s-rrna-gene-sequen-j8nlkpd1g5r7/v2), and for each reverse primer, samples with different forward barcoded primers are pooled together. This script will use regular expressions to locate the forward primer 515F in either R1 or R2 reads, and then split samples based on the barcodes in front of the forward primer.


### requirement
* Python >= 3.5
* logging
* gzip
* argparse


### usage
```
usage: split_samples_by_barcodes_in_fwd_reads.py [-h] [-s START_POSITION]
                                                 [-r PRIMER_REGEX]
                                                 [-b BARCODE_LENGTH]
                                                 [-o OUTPUT_DIR] [-p PREFIX]
                                                 [-l {critical,error,warning,info,debug}]
                                                 [-v]
                                                 barcode_file input_fwd_fastq
                                                 input_rev_fastq

split samples in fastq file based on barcodes found in forward reads at
certain position

positional arguments:
  barcode_file          input barcode file
  input_fwd_fastq       input forward reads in fastq format
  input_rev_fastq       input reverse reads in fastq format

optional arguments:
  -h, --help            show this help message and exit
  -s START_POSITION, --start_position START_POSITION
                        start position of barcodes in sequences, 1-based
                        coordinate, default=5
  -r PRIMER_REGEX, --primer_regex PRIMER_REGEX
                        primer regex, default=GTG[TC]CAGC[AC]GCCGCGGTAA
  -b BARCODE_LENGTH, --barcode_length BARCODE_LENGTH
                        barcode length, default=5
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        output directory, default=./
  -p PREFIX, --prefix PREFIX
                        output prefix, default=split_
  -l {critical,error,warning,info,debug}, --loglevel {critical,error,warning,info,debug}
                        logging level
  -v, --version         show program's version number and exit
```


### example

```bash
#-------------
# Barcode File
#-------------
shengwei@pacific:~/GitHub/515Y926R_Demul/test$ cat A1_barcodes.tsv
#SampleID       BarcodeSequence
A1_1    TCAGC
A1_2    GTATC
A1_3    GCTAC
A1_4    ACGCA
A1_5    GAGAC
A1_6    GACTC
A1_7    CTAGC
A1_8    CGCTC
A1_9    AACAT
```

### Test run

```bash
cd test && ./run_test.sh
```

### count statistics

```bash
#--------
# Results
#--------
shengwei@pacific:~/GitHub/515Y926R_Demul/test/01_demultiplexed_fastq_files$ cat split_barcode2count.tsv 
Barcode	Count
CGCTC	737
CTAGC	643
AACAT	594
GAGAC	579
GACTC	533
GTATC	530
ACGCA	420
GCTAC	345
TCAGC	311
ATCAC	95
NNCTC	12
NNAGC	10
NAGAC	9
NCTAC	8
NNATC	7

shengwei@pacific:~/GitHub/515Y926R_Demul/test/01_demultiplexed_fastq_files$ cat split_sample2count.tsv 
Sample	Count
A1_1	311
A1_2	530
A1_3	345
A1_4	420
A1_5	579
A1_6	533
A1_7	643
A1_8	737
A1_9	594
unknown	308
```
