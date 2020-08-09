#!/bin/bash


# -----------
# demultiplex
# -----------

Demul_Dir="01_demultiplexed_fastq_files"
mkdir -p ${Demul_Dir}

fwd_read=test_rev7.R1.fastq.gz
rev_read=test_rev7.R2.fastq.gz
barcode_file=test_rev7_barcodes.tsv
start_position=5
barcode_length=5
primer_regex="GTG[TC]CAGC[AC]GCCGCGGTAA"
prefix="test_rev7"

#source activate py3k
../split_samples_by_barcodes_in_fwd_reads.py ${barcode_file} ${fwd_read} ${rev_read} \
        -o ${Demul_Dir} -p ${prefix}_ -s ${start_position} \
        -r ${primer_regex} -b ${barcode_length} -l "info" 2>&1 | tee -a 01_demultiplex_stdout_stderr.log
#conda deactivate
gzip ${Demul_Dir}/*.fastq
