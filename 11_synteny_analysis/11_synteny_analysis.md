### Synteny analysis of _Coregonus albula_ assembly

- the 40 chromosomes from our assembly (fCorAlb1) were used for this analysis
- the following _Coregnus_ ssp. genome assemblies were downloaded from NCBI:

  - [_Coregonus_ sp. 'balchen', GCA_902810595.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_902810595.1/)
  - [_Coregonus clupeaformis_, GCF_020615455.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_020615455.1/)
  - [_Coregonus artedi_, GCA_039881085](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_039881085.1/)

- 40 chromosomes were extracted from reference assemblies 'balchen' and _C. clupeaformis_ and the chromosome names were shortened to 1-40
- 38 chromosomes were extracted from reference assembly _C. aretedi_ and the chromosome names were shortened to LG01-LG38

Alignments between the genomes:

```console
C3.fa - our fCorAlb1 assembly
B.fa  - _Coregonus_ sp. 'balchen' assembly
JA.fa - _Coregonus clupeaformis_ assembly
JB.fa - _Coregonus artedi_ assembly
```

were computed with minimap2 (v2.28-r1209).

```bash
minimap2 -x asm10 -t16 C3.fa B.fa | samtools sort -@8 --write-index -o C3_B_asm10.bam
minimap2 -x asm10 -t16 C3.fa JA.fa | samtools sort -@8 --write-index -o C3_JA_asm10.bam
minimap2 -x asm10 -t16 C3.fa JB.fa | samtools sort -@8 --write-index -o C3_JB_asm10.bam
```

JupiterPlot (v1.1) was modified so that it can process bam files instead of sam files:
`jupiterplot/makefile`

```console
# changed sam rule line 86-87 from:
%.bed: %.sam
	grep -E -v '^@' $< | awk '{if($$5 >= $(MAPQ)) print}' | perl $(ROOT_DIR)/bin/samToBed.pl > $@

# into:
%.bed: %.sam
	samtools view $< | grep -E -v '^@' | awk '{if($$5 >= $(MAPQ)) print}' | perl $(ROOT_DIR)/bin/samToBed.pl > $@
```

- run jupiterPlot:

```bash
jupiter name=C3_B_asm10_minBS400K_maxGS400K minBundleSize=400000 gScaff=1 maxGap=400000 ng=0 labels=both ref=C3.fa fa=B.fa t=16 ng=0 labels=both sam=C3_B_asm10.bam
jupiter name=C3_JA_asm10_minBS400K_maxGS400K minBundleSize=400000 gScaff=1 maxGap=400000 ng=0 labels=both ref=C3.fa fa=JA.fa t=16 ng=0 labels=both sam=C3_JA_asm10.bam
jupiter name=C3_JB_asm10_minBS400K_maxGS400K minBundleSize=400000 gScaff=1 maxGap=400000 ng=0 labels=both ref=C3.fa fa=JB.fa t=16 ng=0 labels=both sam=C3_JB_asm10.bam
```

The vector graphic Circos plots (SVG) are available here:
![C3_B_asm10_minBS400K_maxGS400K.svg](C3_B_asm10_minBS400K_maxGS400K.svg)
![C3_JA_asm10_minBS400K_maxGS400K.svg](C3_JA_asm10_minBS400K_maxGS400K.svg)
![C3_JB_asm10_minBS400K_maxGS400K.svg](C3_JB_asm10_minBS400K_maxGS400K.svg)
