# 05_PCA_Admixture_selection_with_LD_pruned_sites

## #PCA 1--all 336 cisco samples

```sh
#!/bin/bash -l
#SBATCH -A uppmax2025-2-114
#SBATCH -p node
#SBATCH -n 1
#SBATCH -C mem256GB
#SBATCH -t 10-00:00:0
#SBATCH -J PCA1_1
#SBATCH -e PCA1_1_%A_%a.err
#SBATCH -o PCA1_1_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=qiaoling.deng@outlook.com

#Load modules
module load bioinfo-tools
module load ANGSD/0.933
#module load python/2.7.15
module load PCAngsd/1.11

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
angsd -bam PCA_feg_sto.txt   \
-ref $REFGENOME -fai $REF_INDEXED \
-remove_bads 1 -only_proper_pairs 1 -minMapQ 30 -minQ 20 \
-GL 2  -doMajorMinor 4 -doMaf 1 -doPost 2 -doGlf 2    \
-sites  $OUT_DIR/ANGSD_Cisco_aassembly2_MAF_005_sites.txt_pruned  \
-out  $OUT_DIR/PCA_Cisco_aassembly2_FegSto_MinMAF_005_pruned  \
-nThreads 8

cd $OUT_DIR

#AdmixRun PCAngsd

for k in {2..10}
do
pcangsd --beagle PCA_Cisco_aassembly2_FegSto_MinMAF_005_pruned.beagle.gz --out ./PCA_Cisco_aassembly2_FegSto_MinMAF_005_pruned  --selection --sites_save  --tree --admix  --admix_K $k  -t 8

done


```

## #PCA2 --Using only samples from Lakes Fegen and Stora Hålsjön

```sh
#!/bin/bash -l
#SBATCH -A uppmax2025-2-114
#SBATCH -p node
#SBATCH -n 1
#SBATCH -C mem256GB
#SBATCH -t 10-00:00:0
#SBATCH -J PCA1_1
#SBATCH -e PCA1_1_%A_%a.err
#SBATCH -o PCA1_1_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=qiaoling.deng@outlook.com

#Load modules
module load bioinfo-tools
module load ANGSD/0.933
#module load python/2.7.15
module load PCAngsd/1.11

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
angsd -bam PCA_feg_sto.txt   \
-ref $REFGENOME -fai $REF_INDEXED \
-remove_bads 1 -only_proper_pairs 1 -minMapQ 30 -minQ 20 \
-GL 2  -doMajorMinor 4 -doMaf 1 -doPost 2 -doGlf 2    \
-sites  $OUT_DIR/ANGSD_Cisco_aassembly2_MAF_005_sites.txt_pruned  \
-out  $OUT_DIR/PCA_Cisco_aassembly2_FegSto_MinMAF_005_pruned  \
-nThreads 8

cd $OUT_DIR

#AdmixRun PCAngsd

for k in {2..10}
do
pcangsd --beagle PCA_Cisco_aassembly2_FegSto_MinMAF_005_pruned.beagle.gz --out ./PCA_Cisco_aassembly2_FegSto_MinMAF_005_pruned  --selection --sites_save  --tree --admix  --admix_K $k  -t 8

done


```

## #PCA3--All other samples after excludingthose from Lakes Fegen and Stora Hålsjön

```sh
#!/bin/bash -l
#SBATCH -A uppmax2025-2-114
#SBATCH -p node
#SBATCH -n 1
#SBATCH -C mem256GB
#SBATCH -t 10-00:00:0
#SBATCH -J PCA2
#SBATCH -e PCA2_%A_%a.err
#SBATCH -o PCA2_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=qiaoling.deng@outlook.com

#Load modules
module load bioinfo-tools
module load ANGSD/0.933
module load PCAngsd/1.11
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
angsd -bam PCA_del_feg_sto.txt  \
-ref $REFGENOME -fai $REF_INDEXED \
-remove_bads 1 -only_proper_pairs 1 -minMapQ 30 -minQ 20 \
-GL 2  -doMajorMinor 4 -doMaf 1 -doPost 2 -doGlf 2    \
-sites  $OUT_DIR/ANGSD_Cisco_aassembly2_MAF_005_sites.txt_pruned  \
-out  $OUT_DIR/PCA_Cisco_aassembly2_del_feg_sto_MinMAF_005_pruned  \
-nThreads 8

cd $OUT_DIR

#AdmixRun PCAngsd

for k in {2..10}
do
pcangsd --beagle PCA_Cisco_aassembly2_del_feg_sto_MinMAF_005_pruned.beagle.gz --out ./PCA_Cisco_aassembly2_del_feg_sto_MinMAF_005_pruned  --selection --sites_save  --tree --admix  --admix_K $k  -t 8
done


```

## #PCA4--Using only Bothnian Bay samples

```sh
#!/bin/bash -l
#SBATCH -A uppmax2025-2-114
#SBATCH -p node
#SBATCH -n 1
#SBATCH -C mem256GB
#SBATCH -t 10-00:00:0
#SBATCH -J PCA3
#SBATCH -e PCA3_%A_%a.err
#SBATCH -o PCA3_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=qiaoling.deng@outlook.com

#Load modules
module load bioinfo-tools
module load ANGSD/0.933
module load PCAngsd/1.11
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
angsd -bam PCA_all_northern.txt  \
-ref $REFGENOME -fai $REF_INDEXED \
-remove_bads 1 -only_proper_pairs 1 -minMapQ 30 -minQ 20 \
-GL 2  -doMajorMinor 4 -doMaf 1 -doPost 2 -doGlf 2    \
-sites  $OUT_DIR/ANGSD_Cisco_aassembly2_MAF_005_sites.txt_pruned  \
-out  $OUT_DIR/PCA_Cisco_aassembly2_all_northern_MinMAF_005_pruned  \
-nThreads 8

cd $OUT_DIR

#AdmixRun PCAngsd

for k in {2..8}
do
pcangsd --beagle PCA_Cisco_aassembly2_all_northern_MinMAF_005_pruned.beagle.gz --out ./PCA_Cisco_aassembly2_all_northern_MinMAF_005_pruned  --selection --sites_save  --tree --admix  --admix_K $k  -t 8

done


```

## #PCA5--Using all Bothnian bay samples after excluding 21 samples from Kalix River  that stand out in panel (d) and two samples from Lule River

```sh
#!/bin/bash -l
#SBATCH -A uppmax2025-2-114
#SBATCH -p node
#SBATCH -n 1
#SBATCH -C mem256GB
#SBATCH -t 10-00:00:0
#SBATCH -J PCA4
#SBATCH -e PCA4_%A_%a.err
#SBATCH -o PCA4_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=qiaoling.deng@outlook.com

#Load modules
module load bioinfo-tools
module load ANGSD/0.933
#module load python/2.7.15
module load PCAngsd/1.11

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
angsd -bam PCA_Brackish_del_LulA19_20.txt  \
-ref $REFGENOME -fai $REF_INDEXED \
-remove_bads 1 -only_proper_pairs 1 -minMapQ 30 -minQ 20 \
-GL 2  -doMajorMinor 4 -doMaf 1 -doPost 2 -doGlf 2    \
-sites  $OUT_DIR/ANGSD_Cisco_aassembly2_MAF_005_sites.txt_pruned  \
-out  $OUT_DIR/PCA_Cisco_aassembly2_Brackish_MinMAF_005_pruned  \
-nThreads 8

cd $OUT_DIR

#AdmixRun PCAngsd

for k in {2..8}
do
pcangsd --beagle PCA_Cisco_aassembly2_Brackish_MinMAF_005_pruned.beagle.gz --out ./PCA_Cisco_aassembly2_Brackish_MinMAF_005_pruned  --selection --sites_save  --tree --admix  --admix_K $k  -t 8

done

```

