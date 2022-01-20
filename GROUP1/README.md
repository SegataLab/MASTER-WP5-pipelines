# MASTER-WP5-Group1 pipeline  

<br>

The workflow is developed by MASTER WP5 group 1 for the analysis of amplicon-based Illumina reads of 16S rRNA gene and ITS sequences. The pipeline processes raw sequencing data and generates amplicon sequence variants (ASVs) gathered in a BIOM-like table, ready for further statistical analysis.  

Authors: [Livio Antonielli](https://github.com/iLivius), [Narciso Martín Quijada](https://github.com/nmquijada).

<br>

## Installation

The workflow runs in a conda environment. Alternatively to [conda](https://github.com/conda/conda), the enviroment can be created with [mamba](https://github.com/mamba-org/mamba). The [pacman](https://github.com/trinker/pacman) package should take care of missing R packages.

```
conda env create -n MASTER-WP5g1 --file DADAism-1.1.0.yml

# or, alternatively:
mamba env create -n MASTER-WP5g1 --file DADAism-1.1.0.yml

conda activate MASTER-WP5g1
```
> Modify ${DADA_PATH} and ${CPUS} variables in `DADAism_bash_master.sh` file.

<br>

## Usage

Launch the `DADAism_bash_master.sh` script within the directory where the FASTQ files are located. Follow the instructions on screen to provide minimal input. In case your amplicons were generated from multiple primers, provide in the same directory a `fwd_primer.fasta` and a `rev_primer.fasta` file, respectively.

<br>

## Workflow description

1) PhiX-contaminant reads, if any, are filtered out with [Bowtie2](https://github.com/BenLangmead/bowtie2). 
2) Primers are stripped off using [Cutadapt](https://github.com/marcelm/cutadapt). 
3) Read length and quality distribution are assessed with [FastQC](https://github.com/s-andrews/FastQC) and inspection of diagnostic plot considered for next steps in [DADA2](https://github.com/benjjneb/dada2).
4) Reads are trimmed and filtered `maxN=0, maxEE=c(2, 2), truncQ=0)`. 
5) Error rate is estimated `randomize=TRUE, MAX_CONSIST=999`, reads are dereplicated and sequence variants inferred. 
6) Denoised forward and reverse reads are merged `minOverlap=10, maxMismatch=1`. 
7) Chimeric ASVs are identified and removed `method="consensus"`
8) Primer-specific V-regions of 16S rRNA gene or ITS regions are targeted and verified using [Metaxa2](https://microbiology.se/software/metaxa2/) or [ITSx](https://microbiology.se/software/itsx/), respectively.
9)  Length distribution of targeted sequences is inspected and refined.
10) Taxonomic assignment of ASVs is carried out with the [RDP Classifier](https://github.com/rdpstaff/classifier) using [SILVA](https://zenodo.org/record/4587955#.YehTM_go_mE) or [UNITE](https://unite.ut.ee/repository.php) as reference `minBoot=80, multithread=TRUE, tryRC=TRUE`.

<br>

## Output

Tables with ASV counts, taxonomic assignment and FASTA files can be found in the `output` directory. Read counts per sample at each processing step and data quality plots can be found in the `stats` directory.