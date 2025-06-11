# 07_1 fst_pi_10kwin_only_signal_freshwater_vs_brackish

## #step1: preparation

```sh
1. canditate gene position
#head List.gene
CA4     RL_17:38000000-42000000
TMTC3   RL_8:51000000-53000000

2. mkdir file for each gene
awk '{print $2 > $1".txt"}' List.gene

3. file name into one file
ls *txt >> canditate_brackish.list
```

## #step2. pi for each group

```sh
#!/bin/bash -l
#SBATCH -A uppmax2025-2-114
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 10-00:00:00
#SBATCH --array=1-7
#SBATCH -J Pi_canditate_brackish
#SBATCH -e Pi_canditate_brackish_%A_%a.err
#SBATCH -o Pi_canditate_brackish_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=qiaoling.deng@outlook.com

#Load modules
module load bioinfo-tools
module load ANGSD/0.933
module load samtools


#STEP 1: Define paths to Refference genome
REFGENOME=/proj/snic2020-2-19/private/cisco/qiaoling/new_assembly/European_cisco_w_extra_data.ipa1.8.0.purged.primary.fa
REF_INDEXED=/proj/snic2020-2-19/private/cisco/qiaoling/new_assembly/European_cisco_w_extra_data.ipa1.8.0.purged.primary.fa.fai


#Go to the bam files directory
BASEDIR=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/Final_Merged_Bam
OUT_DIR=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/Beagle_file/2_MAF_005/0_FST_pi/0_only_signal/3_canditate_gene

#Text file containing sample bam paths
BAM_LIST=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/Beagle_file/2_MAF_005/0_FST_pi

cd $BASEDIR

for POP in `ls $BAM_LIST/Pi_3*.list`; do
POP_TARGET=${POP/.list/}
OUTPUT=$(basename $POP_TARGET)
IND=$(wc -l < "$POP")   #set sample number
MINDEPTH=$(echo "scale=0; $IND * 0.54 / 1" | bc)
MAXDEPTH=$(echo "scale=0; $IND * 1.61 / 1" | bc)

#Region=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/Beagle_file/2_MAF_005/0_FST_pi/0_only_signal/1_spring/LHCGR.txt

GENE_LIST=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/Beagle_file/2_MAF_005/0_FST_pi/0_only_signal/3_canditate_gene/canditate_brackish.list
GENE_NAMES=$(cat $GENE_LIST | sed -n ${SLURM_ARRAY_TASK_ID}p)
GENE=${GENE_NAMES/.txt}
GENE_DIR=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/Beagle_file/2_MAF_005/0_FST_pi/0_only_signal/3_canditate_gene


#Step 3: Finding a 'global estimate' of the SFS per each lake
#------------------------------------------------------------
angsd -bam $POP \
-anc $REFGENOME -fai $REF_INDEXED \
-rf  $GENE_DIR/$GENE_NAMES  \
-doSaf 1 -GL 1 -P 8 -out $OUT_DIR/${OUTPUT}_${GENE}  \
-doCounts 1 -setMinDepth $MINDEPTH  -setMaxDepth $MAXDEPTH  -setMinDepthInd 0.25 -minMapQ 30 -minQ 20 -remove_bads 1 -minInd 10  -skipTriallelic 1 \
-uniqueOnly 1 -dumpCounts 2 -doMaf 1 -doMajorMinor 1
#--------------------------------------------------------------

#Step 3: Obtain the folded site frequency spectrum
realSFS $OUT_DIR/${OUTPUT}_${GENE}.saf.idx -P 8 -fold 1 > $OUT_DIR/${OUTPUT}_${GENE}.sfs

#Step 4: Calculate the thetas for each site
realSFS saf2theta $OUT_DIR/${OUTPUT}_${GENE}.saf.idx  -fold 1 -sfs $OUT_DIR/${OUTPUT}_${GENE}.sfs -outname $OUT_DIR/${OUTPUT}_${GENE}
thetaStat do_stat $OUT_DIR/${OUTPUT}_${GENE}.thetas.idx -win 10000 -step 10000 -outnames $OUT_DIR/${OUTPUT}_${GENE}.theta.10kb.thetasWindow.gz
#Step 5: Estimate Tajimas D and other statistics do a sliding window analysis by adding -win/-step arguments to the last command
#thetaStat print $OUT_DIR/${OUTPUT}_${GENE}.thetas.idx > $OUT_DIR/${OUTPUT}_${GENE}.theta.thetasSNP.gz
#awk '!/^#/ { if ($4 == "-inf") $NF="0"; else $NF=exp($4); print $0 }'  $OUT_DIR/${OUTPUT}_${GENE}.theta.thetasSNP.gz > $OUT_DIR/overlap_${OUTPUT}_${GENE}.theta.thetasSNP.gz

done


```

## #step3: fst (Kalix River--brackish)

```SH
#!/bin/bash -l
#SBATCH -A uppmax2025-2-114
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 10-00:00:00
#SBATCH --array=1-7
#SBATCH -J fst_kalix
#SBATCH -e fst_kalix_%A_%a.err
#SBATCH -o fst_kalix_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=qiaoling.deng@outlook.com

#Load modules
module load bioinfo-tools
module load ANGSD/0.933
module load samtools

## open workdir
List_DIR=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/Pi_Fst/1_list_sh
IN_DIR=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/Beagle_file/2_MAF_005/0_FST_pi/0_only_signal/3_canditate_gene
OUT_DIR=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/Beagle_file/2_MAF_005/0_FST_pi/0_only_signal/3_canditate_gene

FST_LIST_1=$List_DIR/Fst_3_Brackish
FST_LIST_2=$List_DIR/Fst_3_Kalix_21

GENE_LIST=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/Beagle_file/2_MAF_005/0_FST_pi/0_only_signal/3_canditate_gene/canditate_brackish.list
GENE_NAMES=$(cat $GENE_LIST | sed -n ${SLURM_ARRAY_TASK_ID}p)
GENE=${GENE_NAMES/.txt}
GENE_DIR=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/Beagle_file/2_MAF_005/0_FST_pi/0_only_signal/3_canditate_gene


#STEP 2:  calculate the 2dsfs prior
cd $IN_DIR
for group_1 in $(cat $FST_LIST_1); do
  for group_2 in $(cat $FST_LIST_2); do
    echo "Pairs $group_1 and $group_2"
       realSFS $IN_DIR/${group_1}_${GENE}.saf.idx $IN_DIR/${group_2}_${GENE}.saf.idx > $OUT_DIR/${group_1}_${group_2}_${GENE}.ml
    realSFS fst index $IN_DIR/${group_1}_${GENE}.saf.idx $IN_DIR/${group_2}_${GENE}.saf.idx -sfs $OUT_DIR/${group_1}_${group_2}_${GENE}.ml -fstout $OUT_DIR/${group_1}_${group_2}_${GENE}
    realSFS fst stats $OUT_DIR/${group_1}_${group_2}_${GENE}.fst.idx
    realSFS fst stats2 $OUT_DIR/${group_1}_${group_2}_${GENE}.fst.idx -win 10000 -step 10000 > $OUT_DIR/${group_1}_${group_2}_${GENE}_10K.fst_win
    realSFS fst stats2 $OUT_DIR/${group_1}_${group_2}_${GENE}.fst.idx -win 1 -step 1 > $OUT_DIR/${group_1}_${group_2}_${GENE}_fst_nowin
 done
done

```

## #step3: fst (Van+Mal--brackish)

```sh
#!/bin/bash -l
#SBATCH -A uppmax2025-2-114
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 10-00:00:00
#SBATCH --array=1-7
#SBATCH -J fst_brackish
#SBATCH -e fst_brackish_%A_%a.err
#SBATCH -o fst_brackish_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=qiaoling.deng@outlook.com

#Load modules
module load bioinfo-tools
module load ANGSD/0.933
module load samtools

## open workdir
List_DIR=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/Pi_Fst/1_list_sh
IN_DIR=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/Beagle_file/2_MAF_005/0_FST_pi/0_only_signal/3_canditate_gene
OUT_DIR=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/Beagle_file/2_MAF_005/0_FST_pi/0_only_signal/3_canditate_gene

FST_LIST_1=$List_DIR/Fst_4_Brackish
FST_LIST_2=$List_DIR/Fst_4_Van_Sta

GENE_LIST=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/Beagle_file/2_MAF_005/0_FST_pi/0_only_signal/3_canditate_gene/canditate_brackish.list
GENE_NAMES=$(cat $GENE_LIST | sed -n ${SLURM_ARRAY_TASK_ID}p)
GENE=${GENE_NAMES/.txt}
GENE_DIR=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/Beagle_file/2_MAF_005/0_FST_pi/0_only_signal/3_canditate_gene


#STEP 2:  calculate the 2dsfs prior
cd $IN_DIR
for group_1 in $(cat $FST_LIST_1); do
  for group_2 in $(cat $FST_LIST_2); do
    echo "Pairs $group_1 and $group_2"
       realSFS $IN_DIR/${group_1}_${GENE}.saf.idx $IN_DIR/${group_2}_${GENE}.saf.idx > $OUT_DIR/${group_1}_${group_2}_${GENE}.ml
    realSFS fst index $IN_DIR/${group_1}_${GENE}.saf.idx $IN_DIR/${group_2}_${GENE}.saf.idx -sfs $OUT_DIR/${group_1}_${group_2}_${GENE}.ml -fstout $OUT_DIR/${group_1}_${group_2}_${GENE}
    realSFS fst stats $OUT_DIR/${group_1}_${group_2}_${GENE}.fst.idx
    realSFS fst stats2 $OUT_DIR/${group_1}_${group_2}_${GENE}.fst.idx -win 10000 -step 10000 > $OUT_DIR/${group_1}_${group_2}_${GENE}_10K.fst_win
    realSFS fst stats2 $OUT_DIR/${group_1}_${group_2}_${GENE}.fst.idx -win 1 -step 1 > $OUT_DIR/${group_1}_${group_2}_${GENE}_fst_nowin
 done
done


```

