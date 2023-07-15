# Command List

Date: 17-Jul-2023
Author: Moeko Okada

List of commands used in the polyploidy summer school for lectureres.  

## In the morning

## 0. Data download

### Genome Data

```bash
$ cd /tmp/eagle/work/polyploid-rna-seq-analyses/src
src $ wget -O genome.tar.gz https://drive.switch.ch/index.php/s/ZQzcTZ5lGJEcCbA/download
src $ ls
genome.tar.gz  scripts
src $ tar xvzf genome.tar.gz
src $ ls
genome  genome.tar.gz  scripts
src $
```

### downsampled read data

```bash
$ cd /tmp/eagle/work/polyploid-rna-seq-analyses/src
src $ wget -O fastq.tar.gz https://drive.switch.ch/index.php/s/0ny11xweoA5WhEX/download
src $ ls
fastq.tar.gz  genome  genome.tar.gz  scripts
src $ tar xvzf fastq.tar.gz 
src $ ls
fastq  fastq.tar.gz  genome  genome.tar.gz  scripts
src $ ls fastq
MUR_1_R1_ds.fastq.gz  MUR_3_R1_ds.fastq.gz       MUR_R48h_2_R1_ds.fastq.gz  trimmomatic_adapters.txt
MUR_2_R1_ds.fastq.gz  MUR_R48h_1_R1_ds.fastq.gz  MUR_R48h_3_R1_ds.fastq.gz
src $
```

## 1. Quality control

Use shell script to run all samples.  
Make sure you are in the "fastq" directory.

```bash
$ cd /tmp/eagle/work/polyploid-rna-seq-analyses/src/fastq
fastq $ bash ../scripts/1_qc.sh 
fastq $
```

#### 1_qc.sh

```shell
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
```

## 2. Mapping

Use shell script to run all samples.  
Make sure you are in the "1_qc" directory.

```bash
$ cd /tmp/eagle/work/polyploid-rna-seq-analyses/src/1_qc
1_qc $ bash ../scripts/2_map.sh 
1_qc $
```

#### 2_map.sh

```shell
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
```

## 3. Read classification

Use EAGLE-RC software to classify reads.

Use shell script to run all samples.  
Make sure you are in the "2_map" directory.

```bash
$ cd /tmp/eagle/work/polyploid-rna-seq-analyses/src/2_map
2_map $ bash ../scripts/3_eagle.sh 
2_map $
```

#### 3_eagle.sh

```shell
#!/bin/sh
# 3_eagle.sh
# Project: SummerSchool
# Date: 4-Jul-2023
# Author: Moeko Okada

cd /tmp/eagle/work/polyploid-rna-seq-analyses/src
mkdir 3_eagle
cd 2_map/


for i in `ls *_hal.bam`; do
  f=`basename ${i} _hal.bam`
  echo "# start ${f} " $(date)

  eagle-rc --ngi --ref1=../genome/Ahal_v2_2.fa --ref2=../genome/Alyr_v2_2.fa --bam1=../2_map/${i} --bam2=../2_map/${f}_lyr.bam -o ../3_eagle/${f}.eaglerc > ../3_eagle/${f}.eaglerc.stdout.log 2> ../3_eagle/${f}.eaglerc.errout.log

  echo "# done ${f} " $(date)
done
```

## 4. Count reads

Count the number of reads using featureCounts software.

Use shell script to run all samples.  
Make sure you are in the "3_eagle" directory.

```bash
$ cd /tmp/eagle/work/polyploid-rna-seq-analyses/src/3_eagle
3_eagle $ bash ../scripts/4_count.sh 
3_eagle $
```

#### 4_count.sh

```shell
#!/bin/sh
# 4_count.sh
# Project: SummerSchool
# Date: 4-Jul-2023
# Author: Moeko Okada

cd /tmp/eagle/work/polyploid-rna-seq-analyses/src
mkdir 4_count
cd 3_eagle/


for i in `ls *.eaglerc1.ref.bam`; do
  f=`basename ${i} .eaglerc1.ref.bam`
  echo "# start ${f} " $(date)

  featureCounts -T 2 -t exon -g transcript_id -a ../genome/Ahal_v2_2.gtf -o ../4_count/${f}_hal_counts.txt ${i}
  featureCounts -T 2 -t exon -g transcript_id -a ../genome/Alyr_v2_2.gtf -o ../4_count/${f}_lyr_counts.txt ${f}.eaglerc2.ref.bam

  echo "# done ${f} " $(date)
done
```

## 5. Make expression table

```bash
3_eagle $ cd ../4_count
4_count $ pwd
/tmp/eagle/work/polyploid-rna-seq-analyses/src/4_count
4_count $ (echo "gene_id";ls *_hal_counts.txt; ) | sed -e s/_hal_counts.txt//g | tr '\n' '_hal\t' | sed 's/\s*$//' > hal_counts.tsv
4_count $ (echo "gene_id";ls *_lyr_counts.txt; ) | sed -e s/_ds_lyr_counts.txt//g | tr '\n' '_lyr\t' | sed 's/\s*$//' > lyr_counts.tsv
4_count $ echo "" >> hal_counts.tsv
4_count $ echo "" >> lyr_counts.tsv
4_count $ python /tmp/eagle/scripts/tablize.py -skip 1 -a -i 0 -c 6 *_hal_counts.txt >> hal_counts.tsv
4_count $ python /tmp/eagle/scripts/tablize.py -skip 1 -a -i 0 -c 6 *_lyr_counts.txt >> lyr_counts.tsv
4_count $ 
```

## In the afternoon:

## 5. Download the whole expression data

```bash
$ cd src
src $ wget -O exp.tar.gz https://drive.switch.ch/index.php/s/4JaGxRe5nJvJ8WX/download
src $ ls
1_qc  2_map  3_eagle  4_count  exp.tar.gz
fastq  fastq.tar.gz  genome  genome.tar.gz  scripts
src $ tar xvzf exp.tar.gz
src $ ls
1_qc  2_map  3_eagle  4_count  exp  exp.tar.gz
fastq  fastq.tar.gz  genome  genome.tar.gz  scripts
src $ ls exp
hal_counts.tsv      lyr_counts.tsv      metal_genes_ids.txt
src $
```

## 6. Two group comparison using edgeR

```
src $ cd exp
exp $ pwd
/tmp/eagle/work/polyploid-rna-seq-analyses/src/exp
exp $ 

src $ cd 5_deg
5_deg $


```
## 7. Heatmap


## 8. Homoeologous ratio test