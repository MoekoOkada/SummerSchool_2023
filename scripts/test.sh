#!/bin/sh
# test.sh
# Project: SummerSchool
# Date: 4-Jul-2023
# Author: Moeko Okada

PROJECT=/srv/kenlab/moeko/SummerSchool
DS=${PROJECT}/rawdata/downsampled
QC=${PROJECT}/01_qc
EAGLE=${PROJECT}/02_eagle

GENOME=/srv/kenlab/data/genomes/ahal_alyr_genome_v2_2_20160112

eval "$(conda shell.bash hook)"
conda activate ss

# source /usr/local/ngseq/etc/lmod_profile
# module load Tools/samtools/1.17
# module load Tools/EAGLE/1.1.3

fQC=0
if [ ${fQC} -eq 1 ]; then
  cd ${DS}
  echo "# " $(date)

  for i in `ls *_R1_ds.fastq.gz`; do
    f=`basename ${i} _R1_ds.fastq.gz`
    echo "# start ${f} " $(date)

    trimmomatic SE -phred33 -threads 2 \
      ${i} ${QC}/${f}_ds_clean.fastq.gz \
      ILLUMINACLIP:${PROJECT}/trimmomatic_adapters.txt:2:30:10 \
      SLIDINGWINDOW:4:30 MINLEN:40

    echo "# done ${f} " $(date)
  done
fi

fHISAT=0
if [ ${fHISAT} -eq 1 ]; then
  cd ${QC}
  echo "# " $(date)

  for i in `ls *_ds_clean.fastq.gz`; do
    f=`basename ${i} _clean.fastq.gz`
    echo "# start ${f} " $(date)

    hisat2 -p 2 -x ${GENOME}/Ahal_v2_2 -U ${i} | samtools view -SbF4 -q 40 - | samtools sort -m 8G -@ 2 -o ${EAGLE}/${f}_hal.bam
    hisat2 -p 2 -x ${GENOME}/Alyr_v2_2 -U ${i} | samtools view -SbF4 -q 40 - | samtools sort -m 8G -@ 2 -o ${EAGLE}/${f}_lyr.bam
    samtools index ${EAGLE}/${f}_hal.bam
    samtools index ${EAGLE}/${f}_lyr.bam

    echo "# done ${f} " $(date)
  done
fi


fRBH=0
if [ ${fRBH} -eq 0 ]; then
  cd ${GENOME}
  echo "# start genome index " $(date)

  # lastdb -uNEAR -R01 L_origin Alyr_v2_2_1.transcript.CDS.fa
  # lastdb -uNEAR -R01 H_origin Ahal_v2_2_1.transcript.CDS.fa

  lastal L_origin -P2 Ahal_v2_2_1.transcript.CDS.fa | last-map-probs -m 0.49 > L_origin.maf
  lastal H_origin -P2 Alyr_v2_2_1.transcript.CDS.fa | last-map-probs -m 0.49 > H_origin.maf

  gffread -T -o Alyr_v2_2.gtf Alyr_v2_2.gff
  gffread -T -o Ahal_v2_2.gtf Ahal_v2_2.gff

  python ${PROJECT}/homeolog_genotypes.py -f exon -o Ref_L -g Alyr_v2_2.gtf L_origin.maf H_origin.maf
  python ${PROJECT}/homeolog_genotypes.py -f exon -o Ref_H -g Ahal_v2_2.gtf H_origin.maf L_origin.maf

  echo "# finish genome index " $(date)
fi

fCL=0
if [ ${fCL} -eq 0 ]; then
  cd ${EAGLE}
  echo "# " $(date)

  for i in `ls *_lyr.bam`; do
    f=`basename ${i} _lyr.bam`
    echo "# start ${f} " $(date)

    eagle -t 2 -a ${i} -r ${GENOME}/Alyr_v2_2.fa -v ${GENOME}/Ref_L.gtf.vcf --omega=1e-40 --mvh --splice --isc --verbose 1> Ref_L_${f}.txt 2> Ref_L_${f}.readinfo.txt
    eagle-rc -a ${i} --listonly -o Ref_L_${f} Ref_L_${f}.txt Ref_L_${f}.readinfo.txt > Ref_L_${f}.list

    eagle -t 2 -a ${f}_hal.bam -r ${GENOME}/Ahal_v2_2.fa -v ${GENOME}/Ref_H.gtf.vcf --omega=1e-40 --mvh --splice --isc --verbose 1> Ref_H_${f}.txt 2> Ref_H_${f}.readinfo.txt
    eagle-rc -a ${f}_hal.bam --listonly -o Ref_H_${f} Ref_H_${f}.txt Ref_H_${f}.readinfo.txt > Ref_H_${f}.list

    ## Find the consensus classification based on likelihood
    python ${PROJECT}/ref2_consensus.py -A Ref_L_${f}.list -B Ref_H_${f}.list -o ${f}

    ## Write bam files based on consensus list, using A as the reference
    eagle-rc -a ${i} -o ${f}_L --readlist ${f}_L.list
    eagle-rc -a ${f}_hal.bam -o ${f}_H --readlist ${f}_H.list

    ## Perform read counting as you prefer, for example:
    featureCounts -T 2 -t exon -g transcript_id -a ${GENOME}/Alyr_v2_2.gtf -o ${f}_L_counts.txt ${f}_L.bam
    featureCounts -T 2 -t exon -g transcript_id -a ${GENOME}/Ahal_v2_2.gtf -o ${f}_H_counts.txt ${f}_H.bam

    echo "# finish ${f} " $(date)
  done
fi
