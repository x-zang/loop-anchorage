# LoopCore

The LoopCore repository provides a set of related utilities for processing LoopSeq data into assembled contigs. Simple and modular utilites provide

* Simplified installation and deployment
* Easier development of analysis methods

The tools provided in this repository include

* __rust-demultiplex__ : Demultiplex short reads by plate index and sort into bins based on UMI prefix. 
* __rust-umi__ : Given reads binned by UMI prefix, separately group reads by individual UMI sequence. 
* __spades-assembler.py__ : Given a directory of sequence by UMI, use Spades to assembly a contig for each UMI. 
* __main.nf__ : Reference implementation combining the above utilities into end-to-end workflow. 

# Solo

In the various utilties below, some options mention or are specific to the Solo assay. Solo is a modified version of the LoopSeq assay in which each well contains a single long molecule and is indexed with a fixed barcode rather than a random UMI. The structure of this library is also unique, so specific sequence processing steps are required. 

# rust-demultiplex

Demultiplex short reads by plate index and sort into bins based on UMI prefix. 

```
USAGE:
    rust-demultiplex [OPTIONS] <FORWARD_PAIRED> <REVERSE_PAIRED> <FORWARD_UNPAIRED> <BARCODES>

ARGS:
    <FORWARD_PAIRED>      Uncompressed FASTQ file with paired forward reads
    <REVERSE_PAIRED>      Uncompressed FASTQ file with paired reverse reads
    <FORWARD_UNPAIRED>    Uncompressed FASTQ file with unpaired forward reads
    <BARCODES>            comma-delimited list of barcodes for de-multiplexing

OPTIONS:
        --barcode-spacer <BARCODE_SPACER>
            Nucleotide sequence of barcode spacer

    -h, --help
            Print help information

        --pcr-primer <PCR_PRIMER>
            Nucleotide sequence of PCR primer [default: ]

        --pcr-primer-context-length <PCR_PRIMER_CONTEXT_LENGTH>
            Amount of 5' context to use from PCR primer for barcode mapping [default: 2]

        --solo
            Enable processing of Solo data

        --umi-group-length <UMI_GROUP_LENGTH>
            Length of UMI prefix to use for UMI grouping [default: 2]

        --umi-length <UMI_LENGTH>
            Length of the UMI [default: 16]

    -V, --version
            Print version information
```

Additional details

* The rust-demultiplex utility does not perform any trimming of the sequencing adapter. For that reason, sequencing adapter sequences should be trimmed before data is input into rust-demultiplex. The reference main.nf pipeline provides an example of performing such trimming with trimmomatic. 
* UMI group length controls the length of UMI prefix to use for UMI grouping. Longer prefixes will result in more UMI groups. Group length should be selected to roughly match parallel processing available for downstream processing. See nextflow/configs/environments for examples of specific values.
* Barcode spacer sequence, PCR primer sequence and UMI length are fixed arguments based on the assay in question. See nextflow/configs/workflows for examples of specific values. 

# rust-umi

Given reads binned by UMI prefix, separately group reads by individual UMI sequence. 

```
USAGE:
    rust-umi [OPTIONS] <FORWARD_PAIRED> <REVERSE_PAIRED> <FORWARD_UNPAIRED> <SAMPLE_INDEX>

ARGS:
    <FORWARD_PAIRED>      Uncompresed FASTQ file with paired forward reads
    <REVERSE_PAIRED>      Uncompresed FASTQ file with paired reverse reads
    <FORWARD_UNPAIRED>    Uncompresed FASTQ file with unpaired forward reads (currently ignored)
    <SAMPLE_INDEX>        

OPTIONS:
        --barcode-spacer <BARCODE_SPACER>
            Nucleotide sequence of barcode spacer

    -h, --help
            Print help information

        --min-count <MIN_COUNT>
            Minimum number of sequences to output per UMI [default: 100]

        --min-length <MIN_LENGTH>
            Minimum length of trimmed sequence to output [default: 30]

    -p, --pcr-primer <PCR_PRIMER>
            Nucleotide sequence of PCR primer [default: ]

        --solo
            Enable processing of Solo data

    -u, --umi-length <UMI_LENGTH>
            Starting position to search for barcode [default: 16]

    -V, --version
            Print version information
```

Additional details

* In addition, to grouping reads by individual UMI sequence, this tool is also responsible for stripping LoopSeq-specific library sequences. The remaining short reads should represent only the sequence originating from the targeted long read. Therefore, the application of a generic short-read assembly method can recover the target long molecule from the output of this tool.
* Barcode spacer sequence, PCR primer sequence and UMI length are fixed arguments based on the assay in question. See nextflow/configs/workflows for examples of specific values. 


# spades-assembler.py

Given a directory of sequence for UMI, use Spades to assembly a contig for each UMI.

```
usage: Assemble LoopSeq UMI fragments with Spades [-h] --pcr-primer PCR_PRIMER --anchor-start ANCHOR_START --anchor-end ANCHOR_END [--threads THREADS] [--assembly-iterations ASSEMBLY_ITERATIONS]
                                                  [--sampling-target SAMPLING_TARGET] [--solo-contig-barcodes SOLO_CONTIG_BARCODES] [--solo-barcodes SOLO_BARCODES]
                                                  [--r1-orientation R1_ORIENTATION]
                                                  input_dir output_prefix

positional arguments:
  input_dir
  output_prefix

optional arguments:
  -h, --help            show this help message and exit
  --pcr-primer PCR_PRIMER
                        Nucleotide sequence of PCR primer
  --anchor-start ANCHOR_START
                        Expected start sequence of long read
  --anchor-end ANCHOR_END
                        Expected end sequence of long read
  --threads THREADS     Number of threads
  --r1-orientation R1_ORIENTATION
                        Orientation of R1 short reads (if known), relative to long read (FORWARD, REVERSE, or UNKNOWN)
```

Additional Solo options

```
  --assembly-iterations ASSEMBLY_ITERATIONS
                        Number of iterations to repeat assembly (for Solo)
  --sampling-target SAMPLING_TARGET
                        Down-sampling target for each assembly iteration (for Solo)
  --solo-contig-barcodes SOLO_CONTIG_BARCODES
                        Comma-delimited list of contig barcodes (for Solo)
  --solo-barcodes SOLO_BARCODES
                        Comma-delimited list of Solo barcodes
```

Additional details

* This script requires the "pysam" python library. In addition, the tools bwa-mem2 and Spades must be available in the path. 
* If the orientation of the R1 read is known, the orientation of the assembled long read will be updated to the forward direction. The direction is inferred from the alignment of the R1 reads to the original assembled long read. 
* If an anchor start and/or anchor end sequence are provided, these sequences will be trimmed from the start or end of the long read. The start and end sequences are also used to report the status of the assembly (see output below)
* As noted above, certain options only apply to the analysis of Solo data. In analysis of Solo mode, only "UMI" sequences that match expected Solo barcodes will be considered. The short reads for these UMIs will typically provide very high coverage. Assembly is performed in an iterative fashion, where in each iteration reads are down-sampled to the sampling target and assembly is performed. The final assembly is selected as the most frequent full-length result across all iterations (see output below)
* There are two key outputs of this script. The first is a FASTA files containing all of the assembled reads. The second is a CSV file containing some high level general information about the assembled results. The columns of this CSV file are 
  * UMI - The identity of the assembled read as determined by the UMI sequence
  * Length - The length of the assembled read
  * Read count - The number of reads used in the assembly of the long read
  * Status - One of 
    * Full-Length - Both end and start sequences were identified in the assembled read
    * End_only - Only the end sequence was found in the assembled read
    * Start_only - Only the start sequence was found in the assembled read
    * Undetected - Neither start nor end sequence was found in the assembled read

# main.nf

main.nf is a Nextflow workflow that provides a reference implementation of leveraging all of the above utilities to process long reads from unprocessed LoopSeq short read data. The main steps of the workflow are

* Trimming (with Trimmomtic)
* Demultiplex (with rust-demultiplex)
* Read grouping by UMI (with rust-umi)
* Assembly (with Spades/spades_assembler.py)
* Aggregation of results

The Dockerfile in the root directory of this repository can be used to initialize an environment with all of the necessary dependencies to execute the workflow. Nextflow parameters are used to pass options that control both the execution environment and analysis of the workflow.  These options are separated into different Nextflow config files. The "nextflow/configs/environments" directory contains params and config for execution in different computing environments. For example. "devenv.config" is appropriate for running small development runs on a local PC, while "docker.config" would support execution on highly parallel clusters with docker support. The "nextflow/configs/workflows" defines configs with fixed parameter values for certain assay types. Profiles can then be used to easily mix and max different analysis types in different hardware contexts. Profiles are defined in "nextflow.config".

For example, the following command could be used to combined processing of PCR amplicon data in the development environment.  

```
/opt/nextflow/bin/nextflow run main.nf --R1 <input_R1.fastq.gz> --R2 <input_R2.fastq.gz> --sample_list <list of barcodes for demux> --adapters_file <Trimmomatic adapters file> --out_dir <output dir> -profile devenv,pcr_amplicon -c nextflow.config
```

# Development 

This repository makes use of vscode "devcontainers" (https://code.visualstudio.com/docs/devcontainers/containers). With devcontainers configured within vscode, the code in the .devcontainer directory contains all the information needed to initialize a complete development environment for building and running the code in this repository.  