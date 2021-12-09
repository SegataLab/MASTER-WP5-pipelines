# MASTER-WP5-Group6 pipeline

>Authors (alphabetic order): Francesco Rubino (QUB), Narciso Martin Quijada (FFoQSI)

>Version 0.1: initial revision

# Caveats
The current pipeline is a work in progress. There are some aspects that need finalising, including the version of software used. In particular the commands of MGKit pertains to the version 0.5.8, which can be installed using PyPI or conda (using the `frubino` channel). Moreover, in a future version part of these commands will be substituted with Tormes for the functional annotation and later commands will be adapted.

The SNPs calling and pN/pS calculations apply here to a single sample, but in reality needs to be performed on multiple samples at once and that depends on the nature of the experiment. There are portions that can be parallelised, while others need to be applied on multiple samples.

The lines starting with a $ sign (monospace font)  are commands to be executed in a bash shell. The MAGs are assumed to be named `MAG-XX.fa` and the raw reads `SAMPLE-XX-READS-R1.fq.gz` and `SAMPLE-XX-READS-R2.fq.gz`.

Another important note is the naming of the sequences inside each `MAG-XX.fa`, if they come from multiple assemblies. When combining the assemblies it is important to consider a scheme to make the sequence headers unique, for subsequent operations, to integrate data from multiple sources.

# Installation
```bash
conda create -n MASTER-WP5g6 prodigal bowtie2 samtools eggnog-mapper subread pip
pip install --prefix $CONDA_PREFIX git+https://github.com/frubino/mgkit@0.5.8
```
These commands do not include the setup of databases at the moment or further requirements to adapt to specific HPCs.

# Combine MAGs
`cat MAG-*.fa > combined.fa`

# Gene Calling
`prodigal -c -a combined.faa -d combined.fan -f gff -o combined.gff -i combined.fa -p single`

# Using Tormes
There is an option to use TORMES, using a similar approach WP5g7.

```bash
conda deactivate
conda activate WP5g7
echo -e “combinedMAGs\tGENOME\t$(realpath combined.fa) > metadata-tormes.txt
sed -i “1iSamples\tRead1\tRead2” metadata-tormes.txt
tormes -m metadata-tormes.txt -o tormes-combined-MAGs -t ${CPU} \
--only_gene_prediction --prodigal_options “-p single” \
--no_mlst --no_pangenome
conda deactivate
conda activate WP5g6
```
> Substitute `${CPU}` with the number of CPUs you want to use

The `combined.faa` and `combined.gff` that you need will be in `tormes-combined-MAGs/gene_prediction/`

# Functional Annotation
```bash
emapper.py -i combined.faa --output combined -m diamond
Integrate into GFF
edit-gff table -p -ai 20 -c '#' -a FC -t combined.emapper.annotations combined.gff | edit-gff table -p -ai 21 -c '#' -a emapper_desc -t combined.emapper.annotations | edit-gff table -p -ai 1 -c '#' -a emapper_seed -t combined.emapper.annotations | edit-gff table -p -ai 7 -c '#' -a EC -t combined.emapper.annotations | edit-gff table --strip-kegg -p -ai 8 -c '#' -a map_KO -t combined.emapper.annotations | edit-gff table -p -ai 15 -c '#' -a map_CAZY -t combined.emapper.annotations | edit-gff remove -a _dup_score - combined-emapper.gff
```

# Alignment
Assumes that multiple `MAG-XX.fa` files were combined into a single file called `combined.fa` and indexed with `bowtie2-build`:
`bowtie2-build combined.fa combined`

## This is the alignment for a single sample:
```bash
bowtie2 -k 100 --sensitive-local -N 1 --no-unal -p 16 -x combined -1 SAMPLE-XX-READS-R1.fq.gz -2 SAMPLE-XX-READS-R2.fq.gz 2> SAMPLE-XX-ALG.log | samtools view -Sb -F 260 | samtools sort -O bam -T tmp-SAMPLE-XX -o SAMPLE-XX.bam - ; samtools index SAMPLE-XX.bam
```

# Count table
`featureCounts -p -O -f -o combined-counts.tsv -F GTF -a combined.gff -g uid -t CDS`

# Coverage
`samtools depth SAMPLE-XX.bam | alg-utils depth --progress -g combined.gff - SAMPLE-XX.cov`
