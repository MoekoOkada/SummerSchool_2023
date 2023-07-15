# Command List

Date: 17-Jul-2023  
Author: Moeko Okada  

List of commands used in the polyploidy summer school.  

## Part 1

1. Data download
2. Quality control
3. Homoeologous Mapping
4. Read sorting by EAGLE-RC
5. Read count
   1. read count
   2. make a count matrix

## 0. Data download

### Genome Data

```bash
$ cd src
src $ wget -O genome.tar.gz https://drive.switch.ch/index.php/s/ZQzcTZ5lGJEcCbA/download
src $ ls
genome.tar.gz  scripts
src $ tar xvzf genome.tar.gz
src $ ls
genome  genome.tar.gz  scripts
src $
```

### Downsampled read data

```bash
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

Use a shell script to run all samples.  
Make sure you are in the "src" directory.

```bash
src $ bash scripts/1_qc.sh 
src $
```

#### 1_qc.sh: [https://github.com/MoekoOkada/SummerSchool_2023/blob/main/scripts/1_qc.sh](https://github.com/MoekoOkada/SummerSchool_2023/blob/main/scripts/1_qc.sh)

## 2. Mapping

Use a shell script to run all samples.  
Make sure you are in the "src" directory.

```bash
src $ bash scripts/2_map.sh 
src $
```

#### 2_map.sh: [https://github.com/MoekoOkada/SummerSchool_2023/blob/main/scripts/2_map.sh](https://github.com/MoekoOkada/SummerSchool_2023/blob/main/scripts/2_map.sh)


## 3. Read classification

Use EAGLE-RC software to classify reads.

Use a shell script to run all samples.  
Make sure you are in the "src" directory.

```bash
src $ bash scripts/3_eagle.sh 
src $
```

#### 3_eagle.sh: [https://github.com/MoekoOkada/SummerSchool_2023/blob/main/scripts/3_eagle.sh](https://github.com/MoekoOkada/SummerSchool_2023/blob/main/scripts/3_eagle.sh)


## 4. Count reads

Count the number of reads using featureCounts software.

Use a shell script to run all samples.  
Make sure you are in the "src" directory.

```bash
src $ bash scripts/4_count.sh 
src $
```

#### 4_count.sh: [https://github.com/MoekoOkada/SummerSchool_2023/blob/main/scripts/4_count.sh](https://github.com/MoekoOkada/SummerSchool_2023/blob/main/scripts/4_count.sh)


## 5. Make expression table

```bash
src $ cd 4_count
4_count $ (echo "gene_id";ls *_hal_counts.txt; ) | sed -e s/_hal_counts.txt//g | tr '\n' '_hal\t' | sed 's/\s*$//' > hal_counts.tsv
4_count $ (echo "gene_id";ls *_lyr_counts.txt; ) | sed -e s/_ds_lyr_counts.txt//g | tr '\n' '_lyr\t' | sed 's/\s*$//' > lyr_counts.tsv
4_count $ echo "" >> hal_counts.tsv
4_count $ echo "" >> lyr_counts.tsv
4_count $ python /tmp/eagle/scripts/tablize.py -skip 1 -a -i 0 -c 6 *_hal_counts.txt >> hal_counts.tsv
4_count $ python /tmp/eagle/scripts/tablize.py -skip 1 -a -i 0 -c 6 *_lyr_counts.txt >> lyr_counts.tsv
4_count $ 
```
