#!/bin/sh
# 3_eagle.sh
# Project: SummerSchool
# Date: 4-Jul-2023
# Author: Moeko Okada

mkdir 3_eagle
cd 2_map/


for i in `ls *_hal.bam`; do
  f=`basename ${i} _hal.bam`
  echo "# start ${f} " $(date)

  eagle-rc --ngi --ref1=../genome/Ahal_v2_2.fa --ref2=../genome/Alyr_v2_2.fa --bam1=../2_map/${i} --bam2=../2_map/${f}_lyr.bam -o ../3_eagle/${f}.eaglerc > ../3_eagle/${f}.eaglerc.stdout.log 2> ../3_eagle/${f}.eaglerc.errout.log

  echo "# done ${f} " $(date)
done
