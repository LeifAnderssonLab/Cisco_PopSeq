#!/bin/bash -l
#SBATCH -A naiss2023-5-222
#SBATCH -p devcore -n 1
#SBATCH -t 0-00:20:00
#SBATCH -J Cisco_QC
#SBATCH -e /proj/snic2020-2-19/private/cisco/qiaoling/European_Cisco_Genomic_Works/TEMP_outdir_running_logs/cisco_data-%j.err
#SBATCH -o /proj/snic2020-2-19/private/cisco/qiaoling/European_Cisco_Genomic_Works/TEMP_outdir_running_logs/cisco_data-%j.out

#  This script is based on code by Jake Goodall.
# The following script aims to develop the low-pass EU Cisco data from fasta files through to genotype likelihood files.
# The sequencing was done in multiple batches so repeated samples will be combined post mapping

################################################################
##   PHASE 0: Load the appropriate programs for running jobs  ##
################################################################

# Load the programs required for the script
ml load bioinfo-tools
module load bwa/0.7.15
module load ANGSD/0.933
module load samtools/1.12
module load bamtools/2.5.1
module load FastQC/0.11.8
module load MultiQC/1.11
module load QualiMap/2.2.1
module load PCAngsd/0.982
module load trimmomatic/0.39
module load bowtie/1.2.3
module load bowtie2/2.4.5
module load picard/2.23.4
module load GATK/3.8-0
module load GATK/4.2.0.0

################################################################
##  PHASE 1: Index the genome for mapping and setup symlinks  ##
################################################################

### Sample and genome locations for housekeeping

# Individual RAW read data
# /proj/snic2020-2-19/private/cisco/reads/data_from_snpseq00071/snpseq00071/files/VB-3202/221007_A00181_0579_AHYMYLDSX3/

# RUN1 Symlinked Fasta Files
# /proj/snic2020-2-19/private/herring/users/jake/European_Cisco_Library_Quality_and_Preliminary_Structure/RUN1_Symlinks_to_individual_reads

# RUN2 Symlinked Fasta Files
# /proj/snic2020-2-19/private/herring/users/jake/European_Cisco_Library_Quality_and_Preliminary_Structure/RUN2_Symlinks_to_individual_reads

# RUN3 Symlinked Fasta Files
# /proj/snic2020-2-19/private/herring/users/jake/European_Cisco_Library_Quality_and_Preliminary_Structure/RUN3_Symlinks_to_individual_reads

# Cisco assembly location
#/proj/snic2020-2-19/private/cisco/qiaoling/new_assembly



### Generate an index for the cisco genome ###
#cd /proj/snic2020-2-19/private/cisco/qiaoling/new_assembly
#bwa index -p European_cisco_w_extra_data.ipa1.8.0.purged.primary "$cisco_genome"
#bowtie2-build $cisco_genome European_cisco_w_extra_data.ipa1.8.0.purged.primary
#

####################################################
##  PHASE 2: Prepare the fasta files for mapping  ##
####################################################

### Check whether the sample lists are being called properly

# Define variables for easier coding  ###
cisco_genome_bwa_index='/proj/snic2020-2-19/private/cisco/qiaoling/new_assembly/European_cisco_w_extra_data.ipa1.8.0.purged.primary'
cisco_genome_bowtie_index='/proj/snic2020-2-19/private/cisco/qiaoling/new_assembly/European_cisco_w_extra_data.ipa1.8.0.purged.primary'

BASEDIR='/proj/snic2020-2-19/private/herring/users/jake/European_Cisco_Library_Quality_and_Preliminary_Structure'
ADAPTORS='/proj/snic2020-2-19/private/herring/users/jake/European_Cisco_Genomic_Works/FILES_for_processing_samples/Nextera_Mate_Pair_Adaptors.fasta'

SAMPLELIST='/proj/snic2020-2-19/private/herring/users/jake/European_Cisco_Genomic_Works/FILES_for_processing_samples/RUN1_Sample_List_MOD.txt' # Path to the sample list.
RAWFASTQSUFFIX1=_R1_001.fastq.gz # Suffix to raw fastq files. Use forward reads with paired-end data.


# Test whether the sample list is being input properly
#for SAMPLE in `cat $SAMPLELIST`; do
#
#  echo $SAMPLE
#  zcat $BASEDIR'/RUN1_Symlinks_to_individual_reads/'$SAMPLE$RAWFASTQSUFFIX1 | head -n 8
#  echo ' '
#
#done

###  Quality and adaptor trimming of raw reads  ###
###################################################

RAWFASTQDIR=/proj/snic2020-2-19/private/herring/users/jake/European_Cisco_Library_Quality_and_Preliminary_Structure/RUN1_Symlinks_to_individual_reads/ # Path to raw fastq files.
RAWFASTQSUFFIX1=_R1_001.fastq.gz # Suffix to raw fastq files. Use forward reads with paired-end data.
RAWFASTQSUFFIX2=_R2_001.fastq.gz # Suffix to raw fastq files. Use reverse reads with paired-end data.
ADAPTERS='/proj/snic2020-2-19/private/herring/users/jake/European_Cisco_Genomic_Works/FILES_for_processing_samples/Nextera_Mate_Pair_Adaptors.fasta'  # Path to a list of adapter/index sequences.



## Loop over each sample
for SAMPLE in `cat $SAMPLELIST`; do
       SAMPLE_UNIQ_ID=$SAMPLE #define unique sample IDs
       RAWFASTQ_ID=$RAWFASTQDIR$SAMPLE #data input path
       SAMPLEADAPT='/proj/snic2020-2-19/private/cisco/qiaoling/Fasta_Files/RUN_1_trimmed/'$SAMPLE_UNIQ_ID #data output path

       sbatch -J cisco_trim -e /proj/snic2020-2-19/private/cisco/qiaoling/European_Cisco_Genomic_Works/TEMP_outdir_running_logs/$SAMPLE-%j.err -p core -c 1 -t 0-24:00:00 -A naiss2023-5-222 --wrap="trimmomatic PE -threads 1 -phred33 $RAWFASTQ_ID$RAWFASTQSUFFIX1 $RAWFASTQ_ID$RAWFASTQSUFFIX2 \
       $SAMPLEADAPT'_adapter_clipped_f_paired.fastq.gz' $SAMPLEADAPT'_adapter_clipped_f_unpaired.fastq.gz' $SAMPLEADAPT'_adapter_clipped_r_paired.fastq.gz' $SAMPLEADAPT'_adapter_clipped_r_unpaired.fastq.gz' \
       'ILLUMINACLIP:'$ADAPTERS':2:30:10 MINLEN:40'"
done


######################################################
##  PHASE 3: Map fasta to reference and create BAM  ##
######################################################

### Map to the reference, sort, and quality filter  ###
#######################################################

BASEDIR='/proj/snic2020-2-19/private/herring/users/jake/European_Cisco_Library_Quality_and_Preliminary_Structure'

SAMPLELIST=$BASEDIR'/sample_lists/TEST_MAPPING_SCRIPT/TESTMAP_Sample_list_RUN1.txt' # Path to a list of prefixes of the raw fastq files. It should be a subset of the the 1st column of the sample table.
SAMPLETABLE=$BASEDIR'/sample_lists/TEST_MAPPING_SCRIPT/TESTMAP_fastq_table_RUN1.tsv' # Path to a sample table where the 1st column is the prefix of the raw fastq files.
FASTQDIR='/proj/snic2020-2-19/private/cisco/qiaoling/Fasta_Files/RUN_1_trimmed/' # Path to the directory where fastq file are stored.
FASTQSUFFIX1=_adapter_clipped_f_paired.fastq.gz # Suffix to fastq files. Use forward reads with paired-end data.
FASTQSUFFIX2=_adapter_clipped_r_paired.fastq.gz # Suffix to fastq files. Use reverse reads with paired-end data.
MAPPINGPRESET=very-sensitive # The pre-set option to use for mapping in bowtie2 (very-sensitive for end-to-end (global) mapping [typically used when we have a full genome reference]
REFERENCE=/proj/snic2020-2-19/private/cisco/qiaoling/new_assembly/European_cisco_w_extra_data.ipa1.8.0.purged.primary.fa # Path to reference fasta file and file name
BOWTIE_INDEX='/proj/snic2020-2-19/private/cisco/qiaoling/new_assembly/European_cisco_w_extra_data.ipa1.8.0.purged.primary'
REFNAME=European_cisco.hifiasm0.16.1 # Reference name to add to output files, e.g. gadMor2

## Loop over each sample
for SAMPLEFILE in `cat $SAMPLELIST`; do

    ## Extract relevant values from a table of sample, sequencing, and lane ID (here in columns 4, 3, 2, respectively) for each sequenced library
    SAMPLE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 4`
    POP_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 5`
    SEQ_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 3`
    LANE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 2`
    SAMPLE_UNIQ_ID=${SAMPLE_ID}_${POP_ID}_${LANE_ID}  # When a sample has been sequenced in multiple lanes, we need to be able to identify the files from each run uniquely


    ## Extract data type from the sample table
    DATATYPE=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 6`

    ## The input and output path and file prefix
    SAMPLETOMAP=$FASTQDIR$SAMPLE_UNIQ_ID

    SAMPLEBAM='/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/Mapped_Rd1/'$SAMPLE_UNIQ_ID

    ## Define platform unit (PU), which is the lane number
    PU=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 2`

    ## Define reference base name
    REFBASENAME="${REFERENCE%.*}"

    ## Map reads to the reference
    echo $SAMPLE_UNIQ_ID

    # Map the paired-end reads
      if [ $DATATYPE = pe ]; then
    # We ignore the reads that get orphaned during adapter clipping because that is typically a very small proportion of reads. If a large proportion of reads get orphaned (loose their mate so they become single-end), these can be mapped in a separate step and the resulting bam files merged with the paired-end mapped reads.

sbatch -J cisco_trim -e /proj/snic2020-2-19/private/cisco/qiaoling/European_Cisco_Genomic_Works/Temp_Error_Reports/$SAMPLE_UNIQ_ID-%j.err -p core -c 1 -t 10-00:00:00 -A naiss2023-5-222 --wrap="bowtie2 -q --phred33 --$MAPPINGPRESET -p 1 -I 0 -X 1500 --fr --rg-id $SAMPLE_UNIQ_ID --rg SM:$SAMPLE_ID --rg LB:$SAMPLE_ID --rg PU:$PU --rg PL:ILLUMINA -x $BOWTIE_INDEX -1 $SAMPLETOMAP$FASTQSUFFIX1 -2 $SAMPLETOMAP$FASTQSUFFIX2 -S $SAMPLEBAM'_'$SEQ_ID'_'$DATATYPE'_bt2_'$REFNAME'.sam'


    ## Convert to bam file for storage (including all the mapped reads)
    samtools view -bS -F 4 $SAMPLEBAM'_'$SEQ_ID'_'$DATATYPE'_bt2_'$REFNAME'.sam' > $SAMPLEBAM'_'$SEQ_ID'_'$DATATYPE'_bt2_'$REFNAME'.bam'
    #rm -f $SAMPLEBAM'_'$SEQ_ID'_'$DATATYPE'_bt2_'$REFNAME'.sam'

#    ## Filter the mapped reads (to onky retain reads with high mapping quality)
#    # Filter bam files to remove poorly mapped reads (non-unique mappings and mappings with a quality score < 20) -- do we want the quality score filter??
    samtools view -h -q 20 $SAMPLEBAM'_'$SEQ_ID'_'$DATATYPE'_bt2_'$REFNAME'.bam' | samtools view -buS - | samtools sort -o $SAMPLEBAM'_'$SEQ_ID'_'$DATATYPE'_bt2_'$REFNAME'_minq20_sorted.bam'"
    fi
done


######################################################
##  PHASE 4
######################################################
###  Merge all run BAM files together             ###
#####################################################

for SAMPLE in `cat /proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/Merge_Bams/Bam_Group_Names.txt`; do

sbatch -J cisco_merge -e /proj/snic2020-2-19/private/cisco/qiaoling/European_Cisco_Genomic_Works/Temp_Error_Reports/$SAMPLE-%j.err -p core -c 1 -t 5-00:00:00 -A naiss2023-5-222 --wrap="bamtools merge -list '/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/Merge_Bams/'$SAMPLE'.txt' -out '/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/Final_Merged_Bam/'$SAMPLE'_ALL_pe_bt2_European_cisco.hifiasm0.16.1_minq20_sorted.bam'"

done




######################################################
##  PHASE 5
######################################################
###  Deduplicate and clip overlapping read pairs  ###
#####################################################

BAMLIST=/proj/snic2021-6-194/European_Cisco_LAndersson/Data_Files/Bam_Files/Final_Merged_Bam/All_Bam_List_MOD.txt  # It's the same as the fastq setup so we can reuse this file
REFNAME=European_cisco.hifiasm0.16.1 # Reference name to add to output files


## Loop over each sample
for SAMPLEBAM in `cat $BAMLIST`; do
        SAMPLE_UNIQ_ID=$SAMPLEBAM

    ## Remove duplicates and print dupstat file
    ## THIS ONLY NEEDS TO BE DONE IF YOU MERGED THE PE AND SE READS PREVIOUSLY
    # We used to be able to just specify picard.jar on the CBSU server, but now we need to specify the path and version
#    java -jar $PICARD_ROOT/picard.jar MarkDuplicates I=$BASEDIR'/Bam/'$SAMPLEBAM'_bt2_'$REFNAME'_minq20_sorted.bam' O=$BASEDIR'/Bam/'$SAMPLEBAM'_bt2_'$REFNAME'_minq20_sorted_dedup.bam' M=$BASEDIR'/Bam/'$SAMPLEBAM'_bt2_'$REFNAME'_minq20_sorted_dupstat.txt' VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true

   ## Clip overlapping paired end reads (only necessary for paired-end data, so if you're only running se samples, you can comment this step out)
   conda activate cisco_mapping
  sbatch -J cisco_map -e /proj/snic2020-2-19/private/cisco/qiaoling/European_Cisco_Genomic_Works/Temp_Error_Reports/${SAMPLE_UNIQ_ID}%j.err -p core -c 1 -t 1-00:00:00 -A naiss2023-5-222 --wrap="bam clipOverlap --in '/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/Final_Merged_Bam/'$SAMPLEBAM'_bt2_'$REFNAME'_minq20_sorted.bam' --out '/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/Final_Merged_Bam/'$SAMPLEBAM'_bt2_'$REFNAME'_minq20_sorted_overlapclipped.bam' --stats"

done


######################################################
##  PHASE 6   Indel Realignment
######################################################

#!/bin/bash -l
#SBATCH -A naiss2023-5-222
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 10-00:00:00
#SBATCH --array=1-336:1
#SBATCH -J idel_realignment
#SBATCH -e idel_realignment_%A_%a.err
#SBATCH -o idel_realignment_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=qiaoling.deng@imbim.uu.se
#Load modules
module load bioinfo-tools
module load samtools/1.12
module load bamtools/2.5.1
module load GATK-Queue/3.7

#Path to the directory where you have the bam-files
BASEDIR=/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/Final_Merged_Bam/

cd $BASEDIR
# RUN THIS INSIDE THE SEQDIR ALREADY
# Refference genome
REFGENOME='/proj/snic2020-2-19/private/cisco/qiaoling/new_assembly/European_cisco_w_extra_data.ipa1.8.0.purged.primary.fa'

# List the files and get one file:
BAM_FILE=$(ls *_sorted_overlapclipped.bam | sed -n ${SLURM_ARRAY_TASK_ID}p)
OUTPUT=${BAM_FILE/_sorted_overlapclipped.bam/_sorted_overlapclipped_realigned.bam}


# Run programs
java -Xmx38g -jar /sw/bioinfo/GATK/3.7/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R $REFGENOME \
-I $BAM_FILE \
-o $BASEDIR/'all_samples_for_indel_realigner.intervals'


java -Xmx38g -jar /sw/bioinfo/GATK/3.7/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R $REFGENOME \
-I $BAM_FILE \
-targetIntervals $BASEDIR/'all_samples_for_indel_realigner.intervals' \
-o $OUTPUT

# Index all samples running array idexes
  samtools index $OUTPUT



#############################################################################
##  PHASE 7   ###  Check coverage stats with Qualimap multi-bam QC ###
#############################################################################

unset DISPLAY
## Qualimap statistics for all samples individually

BAMQC_SAMPLES_INPUT='/proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/Final_Merged_Bam/EU_Cisco_Sample_List_Clipped_realigned.txt'

qualimap multi-bamqc -data $BAMQC_SAMPLES_INPUT -r -outdir /proj/snic2020-2-19/private/cisco/qiaoling/Bam_Files/Bam/QC -outformat HTML



#############################################################################
##  PHASE 8   ###  visualization  coverage stats with in R ###
#############################################################################
##### copy qualimap output from HTML, generate  remapping.csv

setwd("E:/16-European_cisco/7_remapping")
library("ggplot2")

##load file
cisco_remapping <- read.csv("remapping.csv",header = TRUE)
head(cisco_remapping)
## plot
re_pop <- ggplot(cisco_remapping, aes(x=population, y=Coverage_mean_remapping,)) + 
  stat_boxplot(geom = 'errorbar',linewidth=0.2,cex=1)+
  geom_boxplot(fill = '#8dd3c7' , color = "black")+
  geom_hline(yintercept = 1.0, linetype = "solid", color = "black")+
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red")+
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))
re_pop

re_insert <- ggplot(cisco_remapping, aes(x=population, y=Insert_size_median_remapping)) + 
  stat_boxplot(geom = 'errorbar',linewidth=0.2,cex=1)+
  geom_boxplot(fill = '#bebada' , color = "black")+
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))
re_insert

re_ind <- ggplot(cisco_remapping, aes(x=Sample, y=Coverage_mean_remapping)) + 
  stat_boxplot(geom = 'errorbar',linewidth=0.2,cex=1)+
  geom_boxplot(fill = '#8dd3c7' , color = "black")+
  geom_hline(yintercept = 1.0, linetype = "solid", color = "black")+
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red")+
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))
re_ind

#######Scatter plot
ggplot(data = cisco_remapping, aes(x = Coverage_mean_remapping, y = Insert_size_median_remapping)) +
  geom_point(shape = 18, color = "#377eb8", size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_smooth(method = "lm", se = FALSE, color = "#e41a1c")+
  labs(x = "Mean Coverage", y = "Insert size median") +
  theme_minimal()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))


####### compare
cisco_remapping_compare <- read.csv("remapping_compare.csv",header = TRUE)
head(cisco_remapping_compare)

Coverage <- ggplot(cisco_remapping_compare, aes(x=Mapping, y= Coverage_mean,fill=Mapping)) + 
  stat_boxplot(geom = 'errorbar',linewidth=0.2,cex=1)+
  geom_boxplot()+
  geom_hline(yintercept = 1.0, linetype = "solid", color = "black")+
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red")+
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))
Coverage

Insert <- ggplot(cisco_remapping_compare, aes(x=Mapping, y=Insert_size_median,fill=Mapping)) + 
  stat_boxplot(geom = 'errorbar',linewidth=0.2,cex=1)+
  geom_boxplot()+
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))
Insert



library(gridExtra)
grid.arrange(Coverage, Insert, ncol = 2)


#######plot col
Col_coverage<- ggplot(data = cisco_remapping_compare,aes(x = Sample,y = Coverage_mean,fill = Mapping))+
  geom_col(position = 'dodge',
           width = 0.5)+
  facet_wrap(~Mapping,ncol = 2)+
  theme(axis.text.x = element_blank())+
  geom_hline(yintercept = 1.0, linetype = "solid", color = "black")+
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red")
Col_coverage


Col_Insert<- ggplot(data = cisco_remapping_compare,aes(x = Sample,y = Insert_size_median,fill = Mapping))+
  geom_col(position = 'dodge',
           width = 0.5)+
  facet_wrap(~Mapping,ncol = 2)+
  theme(axis.text.x = element_blank())+
  geom_hline(yintercept = 281.94, linetype = "solid", color = "black")
  
Col_Insert