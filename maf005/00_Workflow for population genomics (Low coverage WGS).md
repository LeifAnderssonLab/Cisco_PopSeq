# Workflow for population genomics (Low coverage WGS)

## 01_mapping

```sh
To get bam-files They must be aligned to the reference, indexed and sorted, named like "id_pop_blablabla.sorted.bam". 

```

## 02_depth_distribution

```sh
This step is for the threshold in the following analysis 
根据distribution计算每个个体的min和max，在后续分析中 min*IND max*IND 来设置depth filiter
```

## 03_run_initial_analysis_on_whole_dataset  -->04,08,09

```sh

this script will work on all bamfiles and calculate saf, maf & genotype likelihood on the whole dataset. It will output in 02_info folder the list of SNP which passed the MIN_MAF and PERCENT_IND filters & their Major-minor alleles (sites_*)


sbatch 01_scripts/01_saf_maf_gl_all.sh
```

## 04_LD pruning with PLINK --->05 ,06, 07

```

```

## 05_PCA_after_LD filter

```

```

## 06_admix_after_LD filter

```

```

## 07_selection_after_LD filter

```

```

## 08_Average_dxy_and_FST ( use first  four chromosome)

```bash
-setMinDepthInd 0.25
```

## 09_association（four contrast） --->  10

```

```

## 10_signal & candidate_gene_set ---> 11,12,13

```

```

## 11_signal_pi  (10kb win) 

```

```

## 12_signal_FST  (Fst is SNP_based)

```

```

## 13_GO_analysis

```

```

