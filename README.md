# K-mer-Based Automated Polishing of T2T Genomes

This session is part of [**Biodiversity Genomics Academy 2024**](https://thebgacademy.org/)

[![Open in Gitpod](https://gitpod.io/button/open-in-gitpod.svg)](https://gitpod.io/#https://github.com/thebgacademy/kmerbased_polishing)

[YouTube Recording](https://www.youtube.com/@thebiodiversitygenomicsacademy) <-- To be updated after BGA24

## Session Leader(s)
Patrick Grady - University of Connecticut

## Description
Introduction to and usage of the updated T2T Automated Polisher from the CHM13 project for application to other genome assemblies.


## Prerequisites

1. Understanding of linux command line basics

>[!warning] "Please make sure you MEET THE PREREQUISITES and READ THE DESCRIPTION above"
>
>    You will get the most out of this session if you meet the prerequisites above.
>
>    Please also read the description carefully to see if this session is relevant to you.
>    
>    If you don't meet the prerequisites or change your mind based on the description or are no longer available at the session time, please email damon at thebgacademy.org to cancel your slot so that someone else on the waitlist might attend.


## Instructions

# T2T-Automated-Polishing

Fully automatic K-mer based polishing of genome assemblies.

Current version is unpublished. Please cite this paper, Arang Rhie's T2T-Polish Git Repository, and this Git Repository if any of the code shared in this repo is used:

Mc Cartney AM, Shafin K, Alonge M et al. Chasing perfection: validation and polishing strategies for telomere-to-telomere genome assemblies. Nat Methods (2022) doi: https://doi.org/10.1038/s41592-022-01440-3

For further details on exact application to T2T-CHM13, read the corresponding section below.

## Description and Best Practices

Auto-Polisher launches an iterative process that allows for more precise K-mer based polishing than typical consensus-only methods. Meryl and Winnowmap2 identify unique k-mers throughout the assembly, and map reads. These reads are extensively filtered in falconc and Merfin to allow for the best base-level polishes. Once corrections are made, this process repeats to now re-anchor on k-mers that are now present in the assembly from previous correction. Generally, base-level accuracy peaks at three iterations (the program default). Genome assembly accuracy can be assessed post-polish by Merqury, and it is highly recommended to use a hybrid k-mer database filtered for k-mers greater than one to obtain the most accurate Merqury QV. The steps for this via [Meryl and Merqury can be found here](https://github.com/arangrhie/T2T-Polish/tree/master/merqury#2-hybrid), as recommended by the developer, Arang Rhie. Using incomplete Meryl DBs to assess post auto-polisher can lead to inaccurate Merqury QV estimates.


## How to Run (Quick Start)

Allocate a fairly large amount of RAM relative to the size of your read set. The Racon step requires the loading of all reads into memory. For instance, a Revio flow cell (~100Gb) requires approximately 400Gb of RAM on a mammalian genome. This pipeline accepts (**and highly recommends**) Herro-corrected ONT reads (use the ONT setting to account for ONT-specific base biases) produced with either the original [Herro](https://github.com/lbcb-sci/herro) (R9 and R10 reads) or [Dorado correct](https://github.com/nanoporetech/dorado) (only R10 at time of v3 release).

```
Automated polishing of draft genomes, version 3

Subcommands:

automated-polishing_v3.sh fullauto - Runs the complete pipeline, including automated GenomeScope k-coverage analysis
automated-polishing_v3.sh polish - Runs a basic automated polishing run
automated-polishing_v3.sh optimizedpolish - Runs an optimized polishing run, requires manual k-cov peak and fitted lookup table from GenomeScope2
automated-polishing_v3.sh computekcov - Calculate kcov and fitted histogram for Merfin using GenomeScope

For help with each subcommand run:
automated-polishing_v3.sh <subcommand> -h|--help
```

```
Fully Automatic Polisher Usage:
 automated-polishing_v3.sh fullauto (options) -d <draft fasta> -r <reads.gz> -s <sequencer type>

Required Arguments:

-d	Draft Fasta - path to the input FASTA/FASTQ containing draft sequences for polishing
-r	Reads - path to the input reads file, in FASTA/FASTQ format (MUST be gzipped).
-s	Sequencing Type - pb or ont, use pb for HiFi and Illumina, ont for all ONT read types

Optional Arguments:

-o	Out Prefix - prefix of the output files. Default: AutoPolisher
-k	K-mer Size - Meryl database k-mer size, should follow initial input. Default (recommended): k=31
-t	Num Threads - number of threads to use. Default: 32
-i	Iterations - number of polishing iterations to perform. Default: 3
```

## Dependencies 
* [Winnowmap2](https://github.com/marbl/Winnowmap)
* [Falconc, available in pbipa package](https://github.com/bio-nim/pb-falconc/releases)
* [Racon (liftover branch)](https://github.com/pgrady1322/racon)
* [Meryl v1.3](https://github.com/marbl/meryl)
* [Merfin v1.0](https://github.com/arangrhie/merfin)
* [Samtools](https://github.com/samtools/samtools)
* [BCFtools](https://github.com/samtools/bcftools)

## Dependencies for FullAuto and ComputeKCov modes
* [Jellyfish](https://github.com/gmarcais/Jellyfish)
* [GenomeScope2](https://github.com/tbenavi1/genomescope2.0)
