BIO373で使っているデータ
  Timの論文
  Kamchatica zink conc. root/leaf
  Expression pattern by RNA-seq
    hal: accumulates zink in leaves and detoxify
    lyr: accumulates zink in roots, and cannot detoxify -> die
  
  HMA4, MTP4 etc... 
    high expression in hal, but low in lyr.

  Tajima's d 
    hal: low diversity => conserved
    lyr: high diversity


  
@fgcz-kl-003:/scratch/bio373_2022/masa/akam_samples2
$ ls /scratch/bio373_2022/data/akam_samples2.tgz
  1 Million lines (250,000 reads) after EAGLE-RC

  HomeoRoq ratioは見せない方向で？


BIO373: thaliana ref
Fribourg: ref => hal, lyr
  /srv/kenlab/data/genomes/ahal_alyr_genome_v2_2_20160112/
  Ahal_v2_2.fa
  Alyr_v2_2.fa
    Convert gene IDs from thaliana to hal/lyr
      oneway_hit_list_Ahal_v2_2_to_TAIR10_cds.xls
      oneway_hit_list_Alyr_v2_2_to_TAIR10_cds.xls
        to find out which gene id is corresponds to HMA etc.

Data will be stored in SwitchDrive, download using wget
  Or GitHub/GitLab



Presentation: qmdが良い


Advanced
  Alignment to thaliana
  (potentially extract reads from bam)
  normal expression analysis using edgeR

Docker from Masa can be used.