# Manual

This repository is an attempt to standardize most of our working
procedures for modern biological experiments. As a side product we
want to achieve real reproducible research, which is one the biggest
problems in modern science.

As we experienced the horror of hundreds of
\*insert-your-favorite-interpreter-language-here\* scripts floating
around without any documentation, we felt the need of a more
sophisticated approach. Thats why we are porting all our workflows to 
[Snakemake](https://bitbucket.org/snakemake/snakemake/wiki/Home).

*Note: Every workflow is a work-in-progress an should be handeld as such. This repository does not make using command-line programs easier. You need to understand how the command-line programs and snakemake work.*

## Required Software

This software is required to be installed in a current version and available via the PATH-variable. Scripts for installation that we use ourselfs are [available here](https://github.com/GrosseLab/InstallProcedures).
* [snakemake](https://bitbucket.org/johanneskoester/snakemake/wiki/Documentation#markdown-header-installation)
* bwa
* segemehl
* bedtools (current and version 2.22.1 - for reasons)
* samtools
* GATK
* picard
* R (with [rpy2](https://pypi.python.org/pypi/rpy2) support)
* featureCounts
* FastQC
* trimmomatic

## List of workflows

* GATK SNP-Calling
* ... other things that you can do by combining all the rules available

## How it should be used

Every set of snakemake rules works only with a fixed naming schema and config-file. Both are documented here.

### File Naming and Directory Structure

Our naming schema uses directories quite heavily to separate different parameters/programs/approaches. This is a subset of the complete naming schema/directory structure, which should suffice to get the gist:

    my_project
    ├── config.json
    ├── Snakefile
    ├── data
    │   └── reads
    │       ├── filtered
    │       └── raw
    │           ├── reads_R1.fastq.gz
    │           └── reads_R2.fastq.gz
    ├── plots
    │   └── coverage
    │       └── bwa
    │           └── hg19
    │               └── raw
    │                   └── BRCA1-ENST00000468300.pdf
    └── results
        ├── coverage
        │   └── bwa
        │       └── hg19
        │           └── raw
        │               └── sample1_illumina_trusightcancer_20_20.cov
        ├── de
        │   └── featureCounts
        │       └── hg19
        │           └── raw
        │               └── all.counts
        ├── mapping
        │   ├── bwa
        │   │   └── hg19
        │   │       └── raw
        │   │           └── sample1.bam
        │   └── segemehl        
        ├── qa
        │   └── fastqc
        │       └── raw
        │           └── sreads_R1_fastqc.html
        └── variants
            └── gatk
                └── hg19
                    └── raw
                        └── patient1.vcf
                        
### Snakefile and config file

The most important part of the config is the definition of the reads-samples-patient relationship. This is done via 2 hashes:

    "samples": {
        "patient1": [ "sample1", "sample2" ],
        "patient2": [ "sampleA", "sampleB" ]
    },
    "units": {
        "sample1": ["reads_R1.fastq.gz", "reads_R2.fastq.gz" ],
        "sample2": ["more_reads_R1.fastq.gz", "more_reads_R2.fastq.gz" ],
        "sampleA": ["S3_R1.fastq.gz", "S3_R2.fastq.gz" ],
        "sampleB": ["S4_R1.fastq.gz", "S4_R2.fastq.gz" ]
    }
    
Sample config and Snakefile s are available in the repository as workflows.

## Release-Numbering conventions

    v0.1p1
     │ │ └─ Patch level: the results are the same but visualization might have
     │ │        changed, bugs preventing output might get fixed
     │ └─── Minor version level: new entire workflows may be available, the
     │          results might have changed due to bug fixes or exchange of
     │          parameters, software, etc.
     └───── Major version level: restructuring of naming schema, results of
                of exitsing workflows might have changed
