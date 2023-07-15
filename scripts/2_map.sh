#!/bin/sh
# 2_map.sh
# Project: SummerSchool
# Date: 4-Jul-2023
# Author: Moeko Okada

cd /tmp/eagle/work/polyploid-rna-seq-analyses/src
mkdir 2_map
cd 1_qc/


for i in `ls *_ds_clean.fastq.gz`; do
  f=`basename ${i} _clean.fastq.gz`
  echo "# start ${f} " $(date)

  hisat2 -p 2 -x ../genome/Ahal_v2_2 -U ${i} | samtools view -SbF4 -q 40 - | samtools sort -m 8G -@ 2 -o ../2_map/${f}_hal.bam
  hisat2 -p 2 -x ../genome/Alyr_v2_2 -U ${i} | samtools view -SbF4 -q 40 - | samtools sort -m 8G -@ 2 -o ../2_map/${f}_lyr.bam
  samtools index ../2_map/${f}_hal.bam
  samtools index ../2_map/${f}_lyr.bam

  echo "# done ${f} " $(date)
done
