#!/bin/sh
# 4_count.sh
# Project: SummerSchool
# Date: 4-Jul-2023
# Author: Moeko Okada

mkdir 4_count
cd 3_eagle/


for i in `ls *.eaglerc1.ref.bam`; do
  f=`basename ${i} .eaglerc1.ref.bam`
  echo "# start ${f} " $(date)

  featureCounts -T 2 -t exon -g transcript_id -a ../genome/Ahal_v2_2.gtf -o ../4_count/${f}_hal_counts.txt ${i}
  featureCounts -T 2 -t exon -g transcript_id -a ../genome/Ahal_v2_2.gtf -o ../4_count/${f}_lyr_counts.txt ${f}.eaglerc1.alt.bam

  echo "# done ${f} " $(date)
done
