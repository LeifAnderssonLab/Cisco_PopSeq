# 03_run_initial_analysis_on_whole_dataset 

## step1: generate beagle file for each chromosome

```bash
#!/bin/bash -l
#SBATCH -A uppmax2025-2-114
#SBATCH -p node
#SBATCH -n 1
#SBATCH -C mem256GB
#SBATCH --array=1-43
#SBATCH -t 1-00:00:0
#SBATCH -J beagle
#SBATCH -e beagle_%A_%a.err
#SBATCH -o begle_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=qiaoling.deng@outlook.com

#Load modules
module load bioinfo-tools
module load ANGSD/0.933
#module load python/2.7.15

######################################################################################
# CONDUCT GENOTYPE LIKELIHOOD ESTIMATION USING ANGSD
# Will output the following files per chromosome:
#  Beagle (genotype likelihood format)


# Qiaoling Deng, 2025/5/14
######################################################################################

#STEP 1: Specify out directory
IN_DIR=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/Final_Merged_Bam
OUT_DIR=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/Beagle_file/2_MAF_005
LIST_DIR=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/List

#STEP 2: Define paths to Reference genome
REFGENOME=/proj/snic2020-2-19/private/cisco/qiaoling/new_assembly/European_cisco_w_extra_data.ipa1.8.0.purged.primary.fa
REF_INDEXED=/proj/snic2020-2-19/private/cisco/qiaoling/new_assembly/European_cisco_w_extra_data.ipa1.8.0.purged.primary.fa.fai

#STEP 3: #Specify chunk list
CHUNK_LIST=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/Beagle_file/2_MAF_005/Cisco_Chunk.list
CHUNK_NAMES=$(cat $CHUNK_LIST | sed -n ${SLURM_ARRAY_TASK_ID}p)
CHUNK_ID=${CHUNK_NAMES/.list}
CHUNK_DIR=/proj/snic2020-2-19/private/cisco/qiaoling/new_assembly/CHR_chunk

###########################################################################################

cd $IN_DIR

#STEP 4: Run ANGSD
angsd -bam $LIST_DIR/list_cisco_all_bam.list \
-ref $REFGENOME -fai $REF_INDEXED \
-rf $CHUNK_DIR/$CHUNK_NAMES \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 20 -minQ 20 \
-GL 2 -trim 0 -doMajorMinor 4 -doMaf 1 -doPost 2 -doGlf 2 -minMaf 0.05 -SNP_pval 1e-6  -docounts 1 \
-setMinDepth 180 -setMaxDepth 540  \
-out $OUT_DIR/Cisco_aassembly2_all_bam_MinMAF_005_${CHUNK_ID} \
-nThreads 8

```

## step2: merge all chromosome

```sh
#!/bin/bash -l
#SBATCH -A uppmax2025-2-114
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 00-10:00:00
#SBATCH -J merge_beagle
#SBATCH -e merge_beagle_%A_%a.err
#SBATCH -o merge_beagle_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=qiaoling.deng@outlook.com


cd /proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/Beagle_file/2_MAF_005

for f in Cisco_aassembly2_all_bam_MinMAF_005_chr_chunk_*.beagle.gz ; do
    zcat $f | tail -n +2 | gzip >> $f.new
done


for f in Cisco_aassembly2_all_bam_MinMAF_005_chr_chunk*gz.new ; do
    cat $f >> Cisco_ANGSD_MAF0.05_Filtered_AllSamples_Output.beagle.gz
done

zcat Cisco_ANGSD_MAF0.05_Filtered_AllSamples_Output.beagle.gz | cut -f1 | rev | sed 's/_/\t/' | rev > ANGSD_Cisco_aassembly2_MAF_005_sites.txt
#######check
head -n 10  ANGSD_Cisco_aassembly2_MAF_005_sites.txt

########index sites
module load bioinfo-tools
module load ANGSD/0.933
angsd sites index  ANGSD_Cisco_aassembly2_MAF_005_sites.txt

```

