#!/bin/bash


# -----------
# demultiplex
# -----------

Demul_Dir="01_demultiplexed_fastq_files"
mkdir -p ${Demul_Dir}

fwd_read=A1_1.fq.gz
rev_read=A1_2.fq.gz 
barcode_file=A1_barcodes.tsv
start_position=1
barcode_length=5
primer_regex="GTG[TC]CAGC[AC]GCCGCGGTAA"
prefix="split"

source activate py3k
python ../split_samples_by_barcodes_in_paired_reads.py ${barcode_file} ${fwd_read} ${rev_read} \
        -o ${Demul_Dir} -p ${prefix}_ -s ${start_position} \
        -r ${primer_regex} -b ${barcode_length} -l "info" 2>&1 | tee -a 01_demultiplex_stdout_stderr.log
conda deactivate
gzip ${Demul_Dir}/*.fastq
