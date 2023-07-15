#!/bin/sh
# 0_ds.sh
# Project: Summer school
# Date: 4-Jul-2023
# Author: Moeko Okada

PROJECT=/srv/kenlab/moeko/SummerSchool
DATA=${PROJECT}/rawdata/roots
DS=${PROJECT}/rawdata/downsampled

cd ${DATA}

for i in `ls *_R1.fastq.gz`; do
  f=${i%%.fastq.gz}
  echo ${f}

  gunzip ${i}
  seqtk sample -s100 ${f}.fastq 1000000 > ${f}_ds.fastq
  gzip ${f}_ds.fastq
  gzip ${f}.fastq
done

