# Genome annotation

The 12_genome_annotation folder is organised in the following folders and files :

```
.
├── maker_abinitio
│   ├── busco
│   │   ├── short_summary.specific.actinopterygii_odb10.busco_evidence_longest_isoform.json
│   │   └── short_summary.specific.actinopterygii_odb10.busco_evidence_longest_isoform.txt
│   └── parameters
│       ├── maker_bopts.ctl
│       ├── maker_evm.ctl
│       ├── maker_exe.ctl
│       └── maker_opts.ctl
├── maker_evidence
│   ├── busco
│   │   ├── short_summary.specific.actinopterygii_odb10.busco_evidence_longest_isoform.json
│   │   └── short_summary.specific.actinopterygii_odb10.busco_evidence_longest_isoform.txt
│   └── parameters
│       ├── maker_bopts.ctl
│       ├── maker_evm.ctl
│       ├── maker_exe.ctl
│       └── maker_opts.ctl
└── statistics
    ├── annotation_stat_after_manual_curation.txt
    ├── annotation_stat_before_manual_curation.txt
    └── final_busco_actinopterygii_odb10.txt
```

# Description of files and folders 

Following is a description of each files and folders :

The maker_abinitio folder contains the busco folder and parameters folder linked to the maker abinitio run (evidence input + augustus input)

```
├── maker_abinitio

```
The busco folder contains the results of busco on the ouput of the maker abinitio run
```
│   ├── busco
│   │   ├── short_summary.specific.actinopterygii_odb10.busco_evidence_longest_isoform.json
│   │   └── short_summary.specific.actinopterygii_odb10.busco_evidence_longest_isoform.txt
```
The parameters folder contains the config files for the maker abinitio run
```
│   └── parameters
│       ├── maker_bopts.ctl
│       ├── maker_evm.ctl
│       ├── maker_exe.ctl
│       └── maker_opts.ctl
```

The maker_evidence folder contains the busco folder and parameters folder linked to the maker evicence run (evidence input meaning proteins, RNAseq)
```
├── maker_evidence
```

The busco folder contains the results of busco on the ouput of the maker evidence run
```
│   ├── busco
│   │   ├── short_summary.specific.actinopterygii_odb10.busco_evidence_longest_isoform.json
│   │   └── short_summary.specific.actinopterygii_odb10.busco_evidence_longest_isoform.txt
```

The parameters folder contains the config files for the maker evidence run
```
│   └── parameters
│       ├── maker_bopts.ctl
│       ├── maker_evm.ctl
│       ├── maker_exe.ctl
│       └── maker_opts.ctl
```

The statistics folder contains the annotation statistics before the manual curation  (before deleting the genes without any functional annotation) and after manual curation. 

It also contains the final busco file, where busco was ran on the final manually curated gff (which is showing the same results as the busco results from the maker abinitio run)

```
└── statistics
    ├── annotation_stat_after_manual_curation.txt
    ├── annotation_stat_before_manual_curation.txt
    └── final_busco_actinopterygii_odb10.txt
```

# Reproducibility
## List of github repositories used and versions during genome annotation analyses

AGAT repository : agat 0.8.1

GAAS repository : gaas 1.2.0

Nextflow pipeline repository: 1.0

## Tools version

Nextflow (22.10.1) singularity-ce (3.8.0) BUSCO (5.4.6) 

fastp (0.23.2) 

hisat2 (2.1.0) 

stringtie(2.2.1) 

RepeatModeler package (2.0.2a) 

RepeatMasker (4.1.2_p1) 

RepeatRunner 

MAKER package (3.01.02) 

exonerate (2.4.0) 

Blast (2.9.0) 

Bioperl (1.7.2) 

Augustus (3.3.3) 

TRNAscan-se (1.3.1) 

Snap (version 2013_11_29) 

GeneMark-ET (4.3) 

GeneMark(ES Suite version 4.48_3.60_lic) 

Interproscan(5.59-91.0)

Infernal (1.1.2) 


## Databases version :

Uniprot Swiss-Prot database (downloaded on 2022-12; 568363 proteins) 

Rfam version 14.9 


