# 01.LD pruning with PLINK

1. purpose

```sh

To remove highly linked SNPs and retain a set of relatively independent SNPs for downstream analyses such as PCA, Admixture, and selection, etc, thereby reducing redundant information and computational burden.

Method:
LD pruning using PLINK’s --indep command. For example:

WINDOW=100: a sliding window of 100 SNPs;

SNP=5: the window shifts by 5 SNPs at a time;

VIF=2: variance inflation factor, which controls the LD threshold (approximately equivalent to r² > 0.5).

Output:
A list of LD-pruned SNPs (a subset with low correlation), which can be used for downstream analyses such as population structure inference.
```

2.tutorial

```sh
https://github.com/clairemerot/angsd_pipeline/tree/master
```

## step 1  plink_LD pruning for the whole genome  (cost 2.5 days to finish)

```sh
#!/bin/bash -l
#SBATCH -A uppmax2025-2-114
#SBATCH -p node
#SBATCH -n 1
#SBATCH -C mem256GB
#SBATCH -t 10-00:00:0
#SBATCH -J plink_prepar
#SBATCH -e plink_prepar_%A_%a.err
#SBATCH -o plink_prepar_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=qiaoling.deng@outlook.com

#Load modules
module load bioinfo-tools
module load ANGSD/0.933
module load plink

#PLINK v1.90b4.9 64-bit (13 Oct 2017)           www.cog-genomics.org/plink/1.9/


######################################################################################
# CONDUCT GENOTYPE LIKELIHOOD ESTIMATION USING ANGSD
# Will output the following files all chromosome:
#  Beagle (genotype likelihood format)


# Qiaoling Deng, 2025/5/14
######################################################################################

#STEP 1: Specify out directory
IN_DIR=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/Final_Merged_Bam
OUT_DIR=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/Beagle_file/2_MAF_005
LIST_DIR=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/List

#STEP 2: Define paths to Refference genome
REFGENOME=/proj/snic2020-2-19/private/cisco/qiaoling/new_assembly/European_cisco_w_extra_data.ipa1.8.0.purged.primary.fa
REF_INDEXED=/proj/snic2020-2-19/private/cisco/qiaoling/new_assembly/European_cisco_w_extra_data.ipa1.8.0.purged.primary.fa.fai

#STEP 3: #Specify chunk list
#CHUNK_LIST=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/Beagle_file/2_MAF_005/Cisco_Chunk.list
#CHUNK_NAMES=$(cat $CHUNK_LIST | sed -n ${SLURM_ARRAY_TASK_ID}p)
#CHUNK_ID=${CHUNK_NAMES/.list}
#CHUNK_DIR=/proj/snic2020-2-19/private/cisco/qiaoling/new_assembly/CHR_chunk

###########################################################################################

cd $IN_DIR

#STEP 4: Run ANGSD
angsd -bam $LIST_DIR/list_cisco_all_bam.list  \
-ref $REFGENOME -fai $REF_INDEXED  \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1  \
-doPlink 2  -doGeno 2 \
-GL 2  -doMajorMinor 4 -doMaf 1 -doPost 2 -doGlf 2 -docounts 1 \
-minMapQ 30 -minQ 20  -SNP_pval 1e-6   -trim 0   \
-sites $OUT_DIR/ANGSD_Cisco_aassembly2_MAF_005_sites.txt \
-out $OUT_DIR/1_plink/Cisco_aassembly2_all_bam_MinMAF_005_plink \
-nThreads 8


####LD prune by plink
OUT_DIR=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/Beagle_file/2_MAF_005

echo "provide a list of LD-pruned snp"
WINDOW=100
SNP=5
VIF=2

cd $OUT_DIR

plink --tped 1_plink/Cisco_aassembly2_all_bam_MinMAF_005_plink.tped \
--tfam 1_plink/Cisco_aassembly2_all_bam_MinMAF_005_plink.tfam \
--indep $WINDOW $SNP $VIF --allow-extra-chr \
--out 1_plink/Cisco_aassembly2_all_bam_MinMAF_005_plink.pruned  \
--threads 8


```

## step2 get the list of pruned sites

```SH
#you need to run first the Rscript to extract the LD-pruned position
#this R script will format the list and angsd will index it
INPUT_plink=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/Beagle_file/2_MAF_005/1_plink/Cisco_aassembly2_all_bam_MinMAF_005_plink.pruned.prune.in
INPUT_angsd=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/Beagle_file/2_MAF_005/ANGSD_Cisco_aassembly2_MAF_005_sites.txt


module load bioinfo-tools
module load R/4.1.1
module load R_packages/4.1.1
module load ANGSD/0.933

Rscript 01_scripts/Rscripts/make_site_list_pruned.r "$INPUT_plink" "$INPUT_angsd"
angsd sites index "$INPUT_angsd"_pruned


```

### Rscript in step2

```R
#this R script pextract the pruned SNPs to make a list of sites fur subsequent analysis in angsd

argv <- commandArgs(T)
INPUT_plink <- argv[1]
INPUT_angsd <- argv[2]

library(dplyr)

pruned<-read.table(INPUT_plink)
head(pruned)
colnames(pruned)<-"LG_pos"
sites<-read.table(INPUT_angsd)
head(sites)
sites$LG_pos<-paste0(sites[,1],"_", sites[,2])


sites_pruned<-inner_join(sites, pruned)
print(paste("there is a total of ", dim(sites)[1], "sites"))
print(paste("we keep ", dim(sites_pruned)[1], "pruned sites"))

write.table(sites_pruned[,1:2], paste0(INPUT_angsd,"_pruned"), col.names=F, row.names=F, quote=F, sep="\t")


```



## step3 run angsd again, with LD pruned sites

```sh
#!/bin/bash -l
#SBATCH -A uppmax2025-2-114
#SBATCH -p node
#SBATCH -n 1
#SBATCH -C mem256GB
#SBATCH -t 10-00:00:0
#SBATCH -J beagle_prune
#SBATCH -e beagle_prune_%A_%a.err
#SBATCH -o begle_prune_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=qiaoling.deng@outlook.com

#Load modules
module load bioinfo-tools
module load ANGSD/0.933
#module load python/2.7.15

######################################################################################
# CONDUCT GENOTYPE LIKELIHOOD ESTIMATION USING ANGSD
# Will output the following files all chromosome:
#  Beagle (genotype likelihood format)


# Qiaoling Deng, 2025/5/20
######################################################################################

#STEP 1: Specify out directory
IN_DIR=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/Final_Merged_Bam
OUT_DIR=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/Beagle_file/2_MAF_005/1_plink
LIST_DIR=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/List

#STEP 2: Define paths to Reference genome
REFGENOME=/proj/snic2020-2-19/private/cisco/qiaoling/new_assembly/European_cisco_w_extra_data.ipa1.8.0.purged.primary.fa
REF_INDEXED=/proj/snic2020-2-19/private/cisco/qiaoling/new_assembly/European_cisco_w_extra_data.ipa1.8.0.purged.primary.fa.fai

cd $IN_DIR

#STEP 4: Run ANGSD
angsd -bam $LIST_DIR/list_cisco_all_bam.list \
-ref $REFGENOME -fai $REF_INDEXED \
-remove_bads 1 -only_proper_pairs 1 -minMapQ 30 -minQ 20 \
-GL 2  -doMajorMinor 4 -doMaf 1 -doPost 2 -doGlf 2    \
-sites  $OUT_DIR/ANGSD_Cisco_aassembly2_MAF_005_sites.txt_pruned  \
-out  $OUT_DIR/Cisco_aassembly2_all_bam_MinMAF_005_pruned  \
-nThreads 8

```

