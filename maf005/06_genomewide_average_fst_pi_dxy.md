# 06_genomewide_average_fst_pi_dxy

#1. pi

```sh
#!/bin/bash -l
#SBATCH -A uppmax2025-2-114
#SBATCH -p node
#SBATCH -n 1
#SBATCH -C mem256GB
#SBATCH -t 10-00:00:00
#SBATCH --array=11,22,33,35
#SBATCH -J Pi
#SBATCH -e Pi_%A_%a.err
#SBATCH -o Pi_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=qiaoling.deng@outlook.com

#Load modules
module load bioinfo-tools
module load ANGSD/0.933
module load samtools



#STEP 1: Define paths to Refference genome
REFGENOME=/proj/snic2020-2-19/private/cisco/qiaoling/new_assembly/European_cisco_w_extra_data.ipa1.8.0.purged.primary.fa
REF_INDEXED=/proj/snic2020-2-19/private/cisco/qiaoling/new_assembly/European_cisco_w_extra_data.ipa1.8.0.purged.primary.fa.fai

#STEP 2: #Specify chunk list
CHUNK_LIST=/proj/snic2020-2-19/private/cisco/qiaoling/new_assembly/CHR_chunk/Cisco_Chunk.list
CHUNK_NAMES=$(cat $CHUNK_LIST | sed -n ${SLURM_ARRAY_TASK_ID}p)
CHUNK_ID=${CHUNK_NAMES/.list}
CHUNK_DIR=/proj/snic2020-2-19/private/cisco/qiaoling/new_assembly/CHR_chunk

#Go to the bam files directory
BASEDIR=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/Final_Merged_Bam
OUT_DIR=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/Beagle_file/2_MAF_005/0_FST_pi

#Text file containing sample bam paths
BAM_LIST=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/Beagle_file/2_MAF_005/0_FST_pi

cd $BASEDIR

for POP in `ls $BAM_LIST/Pi_1*.list`; do
POP_TARGET=${POP/.list/}
OUTPUT=$(basename $POP_TARGET)
IND=$(wc -l < "$POP")   #set sample number
MINDEPTH=$(echo "scale=0; $IND * 0.54 / 1" | bc)
MAXDEPTH=$(echo "scale=0; $IND * 1.61 / 1" | bc)


#Step 3: Finding a 'global estimate' of the SFS per each lake
#------------------------------------------------------------
angsd -bam $POP \
-anc $REFGENOME -fai $REF_INDEXED \
-rf $CHUNK_DIR/$CHUNK_NAMES \
-doSaf 1 -GL 1 -P 8 -out $OUT_DIR/${OUTPUT}_${CHUNK_ID} \
-doCounts 1 -setMinDepth $MINDEPTH  -setMaxDepth $MAXDEPTH  -setMinDepthInd 0.25 -minMapQ 30 -minQ 20 -remove_bads 1 -minInd 10  -skipTriallelic 1 \
-uniqueOnly 1 -dumpCounts 2 -doMaf 1 -doMajorMinor 1
#--------------------------------------------------------------

#Step 3: Obtain the folded site frequency spectrum
realSFS $OUT_DIR/${OUTPUT}_${CHUNK_ID}.saf.idx -P 8 -fold 1 > $OUT_DIR/${OUTPUT}_${CHUNK_ID}.sfs

#Step 4: Calculate the thetas for each site
realSFS saf2theta $OUT_DIR/${OUTPUT}_${CHUNK_ID}.saf.idx  -fold 1 -sfs $OUT_DIR/${OUTPUT}_${CHUNK_ID}.sfs -outname $OUT_DIR/${OUTPUT}_${CHUNK_ID}

#Step 5: Estimate Tajimas D and other statistics do a sliding window analysis by adding -win/-step arguments to the last command

thetaStat do_stat $OUT_DIR/${OUTPUT}_${CHUNK_ID}.thetas.idx -win 10000 -step 10000 -outnames $OUT_DIR/${OUTPUT}_${CHUNK_ID}.theta.10kb.thetasWindow.gz

done

```

#2.fst

```sh
#!/bin/bash -l
#SBATCH -A uppmax2025-2-114
#SBATCH -p node
#SBATCH -n 1
#SBATCH -C mem256GB
#SBATCH -t 6-00:00:00
#SBATCH --array=11,22,33,35
#SBATCH -J fst_1
#SBATCH -e fst_1_%A_%a.err
#SBATCH -o fst_1_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=qiaoling.deng@outlook.com

#-------------------qiaoling 2025-5-30
#Load modules
module load bioinfo-tools
module load ANGSD/0.933
module load samtools

#SET UP CHUNK list
CHUNK_LIST=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/Beagle_file/2_MAF_005/Cisco_Chunk.list
CHUNK_NAMES=$(cat $CHUNK_LIST | sed -n ${SLURM_ARRAY_TASK_ID}p)
CHUNK_ID=${CHUNK_NAMES/.list}
CHUNK_DIR=/proj/snic2020-2-19/private/cisco/qiaoling/new_assembly/CHR_chunk


## open workdir
List_DIR=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/Pi_Fst/1_list_sh
IN_DIR=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/Beagle_file/2_MAF_005/0_FST_pi
OUT_DIR=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/Beagle_file/2_MAF_005/0_FST_pi

FST_LIST_1=$List_DIR/Fst_1_Feg_sto
FST_LIST_2=$List_DIR/Fst_1_del_feg_sto

#STEP 2:  calculate the 2dsfs prior
cd $IN_DIR
for group_1 in $(cat $FST_LIST_1); do
  for group_2 in $(cat $FST_LIST_2); do
    echo "Pairs $group_1 and $group_2"
    for chunk_name in $CHUNK_ID;do
    realSFS $IN_DIR/${group_1}_${chunk_name}.saf.idx $IN_DIR/${group_2}_${chunk_name}.saf.idx > $OUT_DIR/${group_1}_${group_2}_${chunk_name}.ml
    realSFS fst index $IN_DIR/${group_1}_${chunk_name}.saf.idx $IN_DIR/${group_2}_${chunk_name}.saf.idx -sfs $OUT_DIR/${group_1}_${group_2}_${chunk_name}.ml -fstout $OUT_DIR/${group_1}_${group_2}_${chunk_name}
    realSFS fst stats $OUT_DIR/${group_1}_${group_2}_${chunk_name}.fst.idx
    realSFS fst stats2 $OUT_DIR/${group_1}_${group_2}_${chunk_name}.fst.idx -win 10000 -step 10000 > $OUT_DIR/${group_1}_${group_2}_${chunk_name}_10K.fst_win
   done
 done
done

```

3.cal average pi fst dxy--asfsp

```sh
ASFSP=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/Pi_Fst/4_average_fst/asfsp/asfsp.py
module load bioinfo-tools
module load python3
python3 $ASFSP -I  3_fst_output/Pi_2_Autumn_Pi_2_Spring_chr_chunk_26.ml  -D 120,60  --calc dxy

#sfs object is created!
#Number of entries:  [121, 61]
#Pairwise neuclotide difference is:  0.004586714320994426


python3 asfsp.py -I  3_fst_output/Pi_2_Autumn_Pi_2_Spring_chr_chunk_26.ml  -D 120,60  --calc fst
#sfs object is created!
#Number of entries:  [121, 61]
#Hudson's Fst is:  0.11491119436031262


python3 asfsp.py -I  2_Pi_output/Pi_1_Feg_sto_chr_chunk_26.sfs  -D 180  --calc pi
sfs object is created!
#Number of entries:  [181]
#Pairwise neuclotide difference is:  0.004537505903516244

```

