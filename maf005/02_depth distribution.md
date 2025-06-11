# 02_2 depth distribution

## 1. tutorial

```sh
1. https://www.popgen.dk/angsd/index.php/Depth
2. https://github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md#a-tutorial-for-some-basic-analyses-using-ngstoolsangsd

```

## 2. Note

```sh
#Please note that we are analysing only part of chromosome 11, indicated by -r 11. Additionally, -C 50 reduces the effect of reads with excessive mismatches, while -baq 1 computes base alignment quality as explained here (BAQ) to rule out false SNPs close to INDELS, and -trim 0 means that we are not trimming the ends of reads. With -minMapQ 20 we filter out reads with low mapping quality. Finally, -maxDepth 500 means that all sites with depth equal or greater than this value will be binned together, and -P 4 means that I am using 4 threads.
```

## 3. depth distribution

```sh
#!/bin/bash -l
#SBATCH -A uppmax2025-2-114
#SBATCH -p node
#SBATCH -n 1
#SBATCH -C mem256GB
#SBATCH -t 5-00:00:00
#SBATCH --array=11
#SBATCH -J depth_try
#SBATCH -e depth_try_%A_%a.err
#SBATCH -o depth_try_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=qiaoling.deng@outlook.com

#Load modules
module load bioinfo-tools
module load ANGSD/0.933
module load samtools

#STEP 1: Define paths to Refference genome
REFGENOME=/proj/snic2020-2-19/private/cisco/qiaoling/new_assembly/European_cisco_w_extra_data.ipa1.8.0.purged.primary.fa
REF_INDEXED=/proj/snic2020-2-19/private/cisco/qiaoling/new_assembly/European_cisco_w_extra_data.ipa1.8.0.purged.primary.fa.fai

BASEDIR=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/Final_Merged_Bam

OUTDIR=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/Final_Merged_Bam/1_depth_distribution
BAM_LIST=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/Pi_Fst/4_average_fst

#STEP 2: #Specify chunk list
CHUNK_LIST=/proj/snic2020-2-19/private/cisco/qiaoling/new_assembly/CHR_chunk/Cisco_Chunk.list
CHUNK_NAMES=$(cat $CHUNK_LIST | sed -n ${SLURM_ARRAY_TASK_ID}p)
CHUNK_ID=${CHUNK_NAMES/.list}
CHUNK_DIR=/proj/snic2020-2-19/private/cisco/qiaoling/new_assembly/CHR_chunk


for POP in `ls $BAM_LIST/Pi_1_All.list`; do
POP_TARGET=${POP/.list/}
OUTPUT=$(basename $POP_TARGET)

cd $BASEDIR

#We are now printing the distribution of quality scores and per-site depths (global and per-sample).
angsd -P 8 -b $POP -ref $REFGENOME -out $OUTDIR/${OUTPUT}_${CHUNK_ID}_maxdepth.qc -rf $CHUNK_DIR/$CHUNK_NAMES \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0  -C 50  -baq 1 \
        -doDepth 1 -doCounts 1 -maxDepth 2000 -minMapQ 20

done


#in the output .info file, the second colomn is percentile; 
```

## 4.visulization--R

```R
library(ggplot2)
library(ggsci)
library(ggridges)
library(cowplot)
library(scales)
setwd("E:/16-European_cisco/7_review/1-返修/round_1")
depth_chr1 <- read.table(file="Pi_1_All_chr_chunk_1_maxdepth.qc.depthGlobal")
df_t <- as.data.frame(t(depth_chr1))

df_t$new_column <- 0:2000
tail(df_t)
head(df_t)

average_depth <- sum(df_t$new_column * df_t$V1) / sum(df_t$V1)  # ~360
threshold_min <- 0.5 * average_depth #~180
threshold_max <- 1.5 * average_depth #~540

p <- ggplot(df_t, aes(x = new_column, y = V1)) +
  geom_col(width = 1,color = "black") +
  ylim(0,200000)+
  labs(title = " ", x = " ", y = " ") +
  #theme_minimal()
  theme_classic() +
  theme(
    text = element_text(size = 14, face = "bold")  
  ) +
  geom_vline(xintercept = 360, linetype = "dashed", color = "red") + # average
  geom_vline(xintercept = 180, linetype = "dashed", color = "blue") + # min
  geom_vline(xintercept = 540, linetype = "dashed", color = "blue")  # max
 
  

p 


```

5.depth filiter 

```sh
#Threshold: average depth ± 50%

#Note: This depth distribution is based on all samples (n_all = 336). However, when calculating theta , which are computed within each group, the depth filter should be adjusted based on the number of individuals in that specific group (e.g., n_FegS = 30).

#One solution (proposed by Mats) is to scale the global threshold to the group level using the formula: (threshold / n_all) × n_group

threshold_min=0.5 × average depth
threshold_max=1.5 × average depth

-setMinDepth: set (threshold_min / n_all) × n_group  
–setMaxDepth: set (threshold_max / n_all) × n_group  

```

