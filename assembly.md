## Assembly 
```
#!/bin/bash -l
#SBATCH -p node
#SBATCH -N 1
#SBATCH -n 96
#SBATCH -t 2-00:00:00
#SBATCH -J hifiasm_himem.job
#SBATCH -o hifiasm_himem_%j.out
#SBATCH -e hifiasm_himem_%j.err

echo "Start time: `date +%A'-'%k:%M:%S`"

hifiasm -o European_cisco_w_extra_data.hifiasm0.16.1.default -t 96 \
  data/m64077_210809_145516.hifi_reads.fasta.gz \
  data/m64077_210813_140550.hifi_reads.fasta.gz \
  data/m64077_210929_143014.hifi_reads.fasta.gz \
  data/m64204e_210621_010053.hifi_reads.fasta.gz \
  data/m64204e_210626_021411.hifi_reads.fasta.gz \
  data/m64204e_210627_082841.hifi_reads.fasta.gz \
  data/m64204e_210715_115720.hifi_reads.fasta.gz \
  data/m64204e_220209_102511.hifi_reads.fasta.gz \
  data/m64204e_220214_125335.hifi_reads.fasta.gz

echo "End time: `date +%A'-'%k:%M:%S`"
```

## Purge_dups
```
#!/bin/bash -l
#SBATCH -p node
#SBATCH -N 1
#SBATCH -n 96
#SBATCH -t 2-00:00:00
#SBATCH -J purge_dups.job
#SBATCH -o purge_dups_%j.out
#SBATCH -e purge_dups_%j.err

echo "Start time: `date +%A'-'%k:%M:%S`"

singularity exec purge_dups_v1.2.6.sif minimap2 -d European_cisco_w_extra_data.hifiasm0.16.1.default.bp.p_ctg.mmi \
  European_cisco_w_extra_data.hifiasm0.16.1.default.bp.p_ctg.fasta

singularity purge_dups_v1.2.6.sif minimap2 -x map-hifi -I 48G -t 96 European_cisco_w_extra_data.hifiasm0.16.1.default.bp.p_ctg.mmi \
 data/*.fasta.gz  | gzip -c - > European_cisco_w_extra_data.hifiasm0.16.1.default.paf.gz

singularity exec purge_dups_v1.2.6.sif pbcstat European_cisco_w_extra_data.hifiasm0.16.1.default.paf.gz
singularity exec purge_dups_v1.2.6.sif calcuts PB.stat > cutoffs 2>calcults.log
singularity exec purge_dups_v1.2.6.sif hist_plot.py -c cutoffs PB.stat PB.cov.png

singularity exec purge_dups_v1.2.6.sif split_fa \
  European_cisco_w_extra_data.hifiasm0.16.1.default.bp.p_ctg.fasta > European_cisco_w_extra_data.hifiasm0.16.1.default.bp.p_ctg.fasta.split

singularity exec purge_dups_v1.2.6.sif minimap2 -x asm5 -DP -I 48G -t 96 \
  European_cisco_w_extra_data.hifiasm0.16.1.default.bp.p_ctg.fasta.split \
  European_cisco_w_extra_data.hifiasm0.16.1.default.bp.p_ctg.fasta.split \
  | gzip -c - > European_cisco_w_extra_data.hifiasm0.16.1.default.bp.p_ctg.fasta.split.self.paf.gz

singularity exec purge_dups_v1.2.6.sif purge_dups -2 -T cutoffs -c PB.base.cov \
  European_cisco_w_extra_data.hifiasm0.16.1.default.bp.p_ctg.fasta.split.self.paf.gz > dups.bed 2> purge_dups.log

singularity exec purge_dups_v1.2.6.sif get_seqs -p European_cisco_w_extra_data.hifiasm0.16.1.default.bp.p_ctg \
  -e dups.bed European_cisco_w_extra_data.hifiasm0.16.1.default.bp.p_ctg.fasta

echo "End time: `date +%A'-'%k:%M:%S`"
```
