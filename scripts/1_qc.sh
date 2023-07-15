#!/bin/sh
# 1_qc.sh
# Project: SummerSchool
# Date: 4-Jul-2023
# Author: Moeko Okada

cd /tmp/eagle/work/polyploid-rna-seq-analyses/src
mkdir 1_qc
cd fastq/

for i in `ls *_R1_ds.fastq.gz`; do
  f=`basename ${i} _R1_ds.fastq.gz`
  echo "# start ${f} " $(date)

  trimmomatic SE -phred33 -threads 2 \
    ${i} ../1_qc/${f}_ds_clean.fastq.gz \
    ILLUMINACLIP:trimmomatic_adapters.txt:2:30:10 \
    SLIDINGWINDOW:4:30 MINLEN:40

  echo "# done ${f} " $(date)
done
