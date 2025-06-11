# 02_1 depth Manhattan plot

## 1. calculate the average depth across all 336 samples at each site

```sh
#!/bin/bash -l
#SBATCH -A uppmax2025-2-114
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 10-00:00:0
#SBATCH -J depth_all
#SBATCH -e depth_all%A_%a.err
#SBATCH -o depth_all_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=qiaoling.deng@outlook.com

module load bioinfo-tools
module load samtools

REFGENOME=/proj/snic2020-2-19/private/cisco/qiaoling/new_assembly/European_cisco_w_extra_data.ipa1.8.0.purged.primary.fa
REF_INDEXED=/proj/snic2020-2-19/private/cisco/qiaoling/new_assembly/European_cisco_w_extra_data.ipa1.8.0.purged.primary.fa.fai


#STEP 2: Determine directory
BASEDIR=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/Final_Merged_Bam
OUTDIR=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/Final_Merged_Bam/0_depth
BAMLIST=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/Final_Merged_Bam/0_depth/Try.bam.list

cd $BASEDIR

samtools depth -a  -f  $BAMLIST  --reference $REFGENOME  | awk '{ sum = 0; for (i = 3; i <= NF; i++) sum += $i; avg = sum / (NF - 2); print $1, $2, avg; }'  > $OUTDIR/cisco_allsample_average_depth.txt

#-a :Output all positions (including zero depth)
#--reference: Reference sequence FASTA FILE


```

## 2. compute the average sequence coverage in 50 kb windows 

```bash
#!/bin/bash -l
#SBATCH -A uppmax2025-2-114
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 5-00:00:0
#SBATCH -J depth_all
#SBATCH -e depth_all%A_%a.err
#SBATCH -o depth_all_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=qiaoling.deng@outlook.com


#Step1: split file by chromosome

awk '{print > "depth_chr_"$1".txt"}' cisco_allsample_average_depth.txt

#Step2: 50k bin

for chr_file in depth_chr*.txt; do
  awk '{
    win = int($2 / 50000);                            #calculate win order
    key = $1"\t"win*50000"\t"(win+1)*50000;           #make bin label：chr, start, end
    sum[key] += $3;                                  #sum depth across all samples
    count[key] += 1;                                  #count sites in bin
  }
  END {
    for (k in sum) {
      print k"\t"sum[k]"\t"count[k]"\t"sum[k]/count[k];  #output：chr  start  end  sum  count  mean
    }
  }' "$chr_file" | sort -k1,1 -k2,2n > "sorted_${chr_file}"
done


##combine files
cat sorted_depth_chr*.txt > all_chromosomes_sorted_depth.txt


```

## 3. visualization in R

```bash
library(tidyverse)
library(qqman)
library(ggplot2)
library(ggpubr)
setwd("E:/16-European_cisco/7_review/1-返修/round_1")
my_colnames <- c("CHR","BP","winend","snpnumber","count","P")
#P=depth
cisco_depth_all_chr <- read.table(file="all_chromosomes_sorted_depth_del0_del_chrS_chrnumber_for_plot.txt",header=F,col.names=my_colnames)
cisco_depth_all_chr$CHR <- as.numeric(cisco_depth_all_chr$CHR)

head(cisco_depth_all_chr)
cisco_depth_all_chr$SNP <- paste(cisco_depth_all_chr$CHR, cisco_depth_all_chr$BP, sep = "_")
cisco_depth_all_chr_plot <- cisco_depth_all_chr[, c("SNP", "CHR", "BP", "P")]

head(cisco_depth_all_chr_plot)
library(CMplot)
?CMplot

CMplot(cisco_depth_all_chr_plot, 
       
       plot.type = "m", 
       LOG10 = FALSE, 
       col=c("grey", "#d95f02"),
       ylab = "Read Depth", 
       chr.den.col=NULL, cex = 1, axis.cex = 1.5, axis.lwd = 2,
       mar = c(3, 6, 1, 1), lab.cex=2,
       file="tiff",  dpi=1200,
       file.output=TRUE,
       verbose=TRUE, width=24, height=4)


```

