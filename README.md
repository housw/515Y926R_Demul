# 515Y926R_Demul
simple barcode spliter for 515F/926R amplified amplicon sequencing data


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
cd test && ./run_test.sh
```

### count statistics

```bash
$ cd test/01_demultiplexed_fastq_files
$ head test_rev7_barcode2count.tsv 
Barcode	Count
CGCGT	1181
CGCTA	901
AGCTT	885
TCAGC	853
GCTAC	716
GTATC	684
GGAAC	680
CCCTT	657
ACAAC	654
$ head test_rev7_sample2count.tsv 
Sample	Count
A1A_a	647
A1B_a	493
A1C_a	574
A2A_a	606
A2B_a	660
A2C_a	606
A3A_a	519
T0_a	606
T0_b	908
```
