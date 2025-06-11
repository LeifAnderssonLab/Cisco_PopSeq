# 10_consensus_fasta

```sh
#!/bin/bash -l
#SBATCH -A uppmax2025-2-114
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 00-10:00:00
#SBATCH -J consensus_pop
#SBATCH -e consensus_pop_%A_%a.err
#SBATCH -o consensus_pop_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=qiaoling.deng@outlook.com

#Load modules
module load bioinfo-tools
module load samtools

####### Qiaoling 2025-6-3

#Go to the bam files directory
BASEDIR=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/Final_Merged_Bam
OUT_DIR=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/ANGSD_analysis/Beagle_file/2_MAF_005/4_consensus

cd $BASEDIR

#### merge all 2 south population BAM file
# -b FILE: List of input BAM files, one file per line.

samtools merge -o $OUT_DIR/merged_FegS.bam FegS*_overlapclipped_realigned.bam

#### create bam index
cd $OUT_DIR

samtools index -M $OUT_DIR/merged_*.bam

samtools consensus  -r RL_26:15725105-15728166 -a --show-ins no merged_FegS.bam  -o merged_FegS_consensus_LHCGR.fa


```

