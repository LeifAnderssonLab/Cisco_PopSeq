# 09 maf_pop

## #step1: maf  for each pop

```sh
#!/bin/bash -l
#SBATCH -A uppmax2025-2-114
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 10-00:00:00
#SBATCH --array=1-13
#SBATCH -J maf
#SBATCH -e maf_%A_%a.err
#SBATCH -o maf_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=qiaoling.deng@outlook.com

#Load modules
module load bioinfo-tools
module load ANGSD/0.933

######################################################################################
# Qiaoling Deng,2025-6-3
######################################################################################

#Specify which bam list to use
BAM_LIST_PATH=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/Beagle_file/2_MAF_005/3_maf_pop


#Specify which bam to use
BAM_LIST=$(ls $BAM_LIST_PATH/Maf*.list| sed -n ${SLURM_ARRAY_TASK_ID}p)
BAM_TARGET=${BAM_LIST/.list/}
OUTPUT=$(basename $BAM_TARGET)

#Sites
SITES=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/Beagle_file/2_MAF_005/ANGSD_Cisco_aassembly2_MAF_005_sites.txt

#Specify the number of individuals
x=`cat $BAM_LIST | wc -l`


#STEP 1: Define paths to Reference genome
REFGENOME=/proj/snic2020-2-19/private/cisco/qiaoling/new_assembly/European_cisco_w_extra_data.ipa1.8.0.purged.primary.fa
REF_INDEXED=/proj/snic2020-2-19/private/cisco/qiaoling/new_assembly/European_cisco_w_extra_data.ipa1.8.0.purged.primary.fa.fai


#Out directory
OUT_DIR=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/Beagle_file/2_MAF_005/3_maf_pop

#Go to the bam files directory

BASEDIR=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/Final_Merged_Bam
cd $BASEDIR

#STEP 9: The association can then be performed on the genotype probabilities using the score statistics
angsd -out $OUT_DIR/$OUTPUT \
-ref $REFGENOME -fai $REF_INDEXED \
-doMajorMinor 4 -doMaf 1 -bam $BAM_LIST \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 30 -minQ 20 \
-gl 2 \
-sites $SITES \
-nThreads 8


```

#step 2: combine

```sh
#!/bin/bash -l
#SBATCH -A uppmax2025-2-114
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 10-00:00:00
#SBATCH -J maf_zip
#SBATCH -e maf_zip_%A_%a.err
#SBATCH -o maf_zip_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=qiaoling.deng@outlook.com


#Load modules
module load bioinfo-tools
module load ANGSD/0.933
#step1:解压缩

gunzip Maf_*.gz

#step2：加入群体信息
awk 'BEGIN {OFS="\t"} {print $0, "kalixA"}' Maf_kalix_21.mafs > Maf_pop_kalixA.mafs
awk 'BEGIN {OFS="\t"} {print $0, "kalix_hybrid"}' Maf_north_kalixA_9.mafs > Maf_pop_north_kalix_hybrid.mafs
awk 'BEGIN {OFS="\t"} {print $0, "kalixK"}' Maf_north_kalk.mafs > Maf_pop_north_kalixK.mafs
awk 'BEGIN {OFS="\t"} {print $0, "LulA"}' Maf_north_LulA.mafs > Maf_pop_north_LulA.mafs
awk 'BEGIN {OFS="\t"} {print $0, "PitA"}' Maf_north_PitA.mafs > Maf_pop_north_PitA.mafs
awk 'BEGIN {OFS="\t"} {print $0, "PitK"}' Maf_north_PitK.mafs > Maf_pop_north_PitK.mafs
awk 'BEGIN {OFS="\t"} {print $0, "UleO"}' Maf_north_UleO.mafs > Maf_pop_north_UleO.mafs

awk 'BEGIN {OFS="\t"} {print $0, "FegA"}' Maf_south_FegA.mafs > Maf_pop_south_FegA.mafs
awk 'BEGIN {OFS="\t"} {print $0, "FegS"}' Maf_south_FegS.mafs > Maf_pop_south_FegS.mafs
awk 'BEGIN {OFS="\t"} {print $0, "Sta"}' Maf_south_Sta.mafs > Maf_pop_south_Sta.mafs
awk 'BEGIN {OFS="\t"} {print $0, "Sto"}' Maf_south_Sto.mafs > Maf_pop_south_Sto.mafs
awk 'BEGIN {OFS="\t"} {print $0, "VanO"}' Maf_south_VanO.mafs > Maf_pop_south_VanO.mafs
awk 'BEGIN {OFS="\t"} {print $0, "VanV"}' Maf_south_VanV.mafs > Maf_pop_south_VanV.mafs

#step3：创建一个名为 "Maf_pop_all_cisco.mafs" 的空文件，作为合并后的结果文件
touch Maf_pop_all_cisco.mafs

#step4：合并文件,不含第一行
for file in Maf_pop*.mafs;do
 tail -n +2 $file >> Maf_pop_all_cisco.mafs
done

#step5：压缩并删除中间文件
gzip *.mafs

```

