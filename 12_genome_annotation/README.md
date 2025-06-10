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
