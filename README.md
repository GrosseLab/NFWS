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

## Installation
In an earlier version we used Docker to ensure the consistency between pipeline installtions. Today we decided to use the conda package manager along the bioconda and conda-forge channels.

1. Please first install miniconda on your machine using these resources: https://conda.io/en/latest/miniconda.html.
2. Afterwards, please setup the bioconda and conda-forge channels as described here:https://bioconda.github.io/#set-up-channels.
3. Next, you should install snakemake and some useful helper tools in a new environment
```
conda create -n NFWS snakemake git pigz samtools
conda activate snakemake
```
4. You can now execute the NFWS snakemake pipeline with conda
```
snakemake --use-conda
```

## List of analyses

* GATK4 SNP-Calling
* CNV analysis
* coverage analysis
* Star-salmon-tximport-edgeR RNA-Seq analysis
* ... other things that you can do by combining all the rules available

## How it should be used

Every set of snakemake rules works only with a fixed naming schema and config-file. Both are documented here.

### File Naming and Directory Structure

Our naming schema uses directories quite heavily to separate different parameters/programs/approaches. This is a subset of the complete naming schema/directory structure, which should suffice to get the gist:

    my_project
    ├── config.yaml
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

The most important part of the config is the definition of the reads-samples-patient relationship. This is done via 3 hashes for condition, biological replicates, and technical replicates respectively:

    replicates:
      treatment: [patient1, patient2]
      control: [patient3, patient4]
    samples:
      patient1: [ sample1A, sample1B ]
      patient2: [ sample2A, sample2B ]
      patient3: [ sample3A, sample3B ]
      patient4: [ sample4A, sample4B ]
    },
    units:
      sample1A: [reads_R1.fastq.gz, reads_R2.fastq.gz ]
      sample1B: [more_reads_R1.fastq.gz, more_reads_R2.fastq.gz ]
      sample2A: [S3_R1.fastq.gz, S3_R2.fastq.gz ]
      sample2B: [S4_R1.fastq.gz, S4_R2.fastq.gz ]
      sample2A: [S5_R1.fastq.gz, S5_R2.fastq.gz ]
      sample2B: [S6_R1.fastq.gz, S6_R2.fastq.gz ]
      sample2A: [S7_R1.fastq.gz, S7_R2.fastq.gz ]
      sample2B: [S8_R1.fastq.gz, S8_R2.fastq.gz ] 
    }
    
Sample config and Snakefiles are available in the repository as workflows. Some may still contain json instead of yaml config files.

### Parameter configuration
Parameters for programs can be specified in the config file. They are organized by parameter sets (i.e. code names for a set of parameters for several programs) and program names:

    parameters:
      gtf_adj:
        star: ["--sjdbGTFfeatureExon", "exon", "--sjdbGTFtagExonParentTranscript", "Parent"]
        salmon: ["--fldMean",  "190"]

## Release-Numbering conventions

    v0.1p1
     │ │ └─ Patch level: the results are the same but visualization might have
     │ │        changed, bugs preventing output might get fixed
     │ └─── Minor version level: new entire workflows may be available, the
     │          results might have changed due to bug fixes or exchange of
     │          parameters, software, etc.
     └───── Major version level: restructuring of naming schema, results of
                of exitsing workflows might have changed
