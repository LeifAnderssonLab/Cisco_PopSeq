### Manual curation of _Coregonus albula_ assembly

The manual curation was done in HiC contact maps:

- high resolution pretext files (\*.pretext) which can be visualizes and modified with [PretextView](https://github.com/sanger-tol/PretextView)
- multi-resolution cooler files (\*.mcool) which can visualized with [HiGlass](https://github.com/higlass/higlass-docker)

pretext and mcool files were generated with the [Earth-Biogenome-Project-pilot](https://github.com/NBISweden/Earth-Biogenome-Project-pilot/tags) assembly nextflow pipeline (`tag: Coregonus-albula_curation`).

workflow_parameters.yaml

```yaml
project: "uppmax2025-2-58"
input: "assembly_parameters.yml"
steps: "curation"
telomer_motif: "AACCCT"
organelle_assembly_mode: none
```

assembly_parameters.yaml

```yaml
# Mandatory - sample metadata
sample:
  name: "Coregonus albula"

hifi:
  - reads: "../../data/raw-data/PacBio-WGS/m64077_210809_145516.hifi_reads.fasta"
  - reads: "../../data/raw-data/PacBio-WGS/m64077_210813_140550.hifi_reads.fasta"
  - reads: "../../data/raw-data/PacBio-WGS/m64077_210929_143014.hifi_reads.fasta"
  - reads: "../../data/raw-data/PacBio-WGS/m64204e_210621_010053.hifi_reads.fasta"
  - reads: "../../data/raw-data/PacBio-WGS/m64204e_210626_021411.hifi_reads.fasta"
  - reads: "../../data/raw-data/PacBio-WGS/m64204e_210627_082841.hifi_reads.fasta"
  - reads: "../../data/raw-data/PacBio-WGS/m64204e_210715_115720.hifi_reads.fasta"
  - reads: "../../data/raw-data/PacBio-WGS/m64204e_220209_102511.hifi_reads.fasta"
  - reads: "../../data/raw-data/PacBio-WGS/m64204e_220214_125335.hifi_reads.fasta"

hic:
  - read1: "../../data/raw-data/Illumina-HiC/P21653_201_S6_L001_R1_001.fastq.gz"
    read2: "../../data/raw-data/Illumina-HiC/P21653_201_S6_L001_R2_001.fastq.gz"

assembly:
  - assembler: "hifiasm"
    stage: "scaffolded" # available stages: 'raw', 'decontaminated', 'polished', 'scaffolded', 'curated'
    id: "fCorAlb_scaffolds"
    pri_fasta: "../../data/asm/fCorAlb_scaffolds.fasta"
```

The pipeline was started with:

```bash
nextflow run "NBISweden/Earth-Biogenome-Project-pilot" \
        -r "Coregonus-albula_curation" \
        -profile "rackham" \
        -work-dir "nxf-work" \
        -resume \
        -ansi-log false \
        -params-file workflow_parameters.yml \
        --cache "nobackup/database-cache" \
        --outdir "results"
```

The high resolution pretext map including the `gap-`,`coverage-` and `telomere-track` were opened in PretextView.
The mcool and the contig grid were ingested into a local HiGlass-Docker instance.

```bash
## ingest multi-cooler file into higlass
project=fCorAlb_mc01;
f=fCorAlb.dups.mcool
cp $f ~/hg-tmp && docker exec higlass-container python higlass-server/manage.py ingest_tileset --filename /tmp/${f} --filetype cooler --datatype matrix --project-name ${project} && rm ~/hg-tmp/$f
## ingest contig grid file into higlass
project=fCorAlb_mc01;
f=fCorAlb_final_asm.contigs.sizes
cp $f ~/hg-tmp && docker exec higlass-container python higlass-server/manage.py ingest_tileset --filename /tmp/${f} --filetype chromsizes-tsv --datatype chromsizes --project-name ${project} && rm ~/hg-tmp/$f
```

In total two manual curation rounds were performed and which include 276 interventions per Gb.
