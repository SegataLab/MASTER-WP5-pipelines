# MASTER-WP5-Group7 pipeline  

<br>

This is the pipeline developed by MASTER WP5 group 7 for the screening of antimicrobial resistance (AMR) and virulence genes. Additionally, it will generate the gene prediction files (*gff, faa, fna*) that will be used for the metagenomic annotation by MASTER WP5 group 6.  
The pipeline can be run in two modes depending on the input data: [assembly-based](#running-assembly-based-pipeline) and [assembly-free](#running-assembly-free-pipeline).  

Authors (order as in the "MASTER_BIOINF Tools" GoogleDoc): [Narciso Martín Quijada](https://github.com/nmquijada), [José F. Cobo-Díaz](https://github.com/JoseCoboDiaz), Francesca De Filippis, [Raúl Cabrera-Rubio](https://github.com/RaulCR), [Carlos Sabater](https://github.com/CarlosSabaterSanchez) and [Francesco Rubino](https://github.com/frubino).

<br>

## Installation

The pipeline has been devised to be run as a conda environment.

```
wget https://anaconda.org/nmquijada/tormes-1.3.0/2021.06.08.113021/download/tormes-1.3.0.yml 

conda env create -n MASTER-WP5g7 --file tormes-1.3.0.yml

conda activate MASTER-WP5g7
tormes-setup

# Get bowtie2
conda install bowtie2
bowtie2-build ${PATH_TO_CONDA_ENVS}/MASTER-WP5g7/db/resfinder/sequences \
  ${PATH_TO_CONDA_ENVS}/MASTER-WP5g7/db/resfinder/resfinder.bowtie --threads ${CPUS}
bowtie2-build ${PATH_TO_CONDA_ENVS}/MASTER-WP5g7/db/vfdb/sequences \
  ${PATH_TO_CONDA_ENVS}/MASTER-WP5g7/db/vfdb/vfdb.bowtie --threads ${CPUS}
bowtie2-build ${PATH_TO_CONDA_ENVS}/MASTER-WP5g7/db/card/sequences \
  ${PATH_TO_CONDA_ENVS}/MASTER-WP5g7/db/card/card.bowtie --threads ${CPUS}

# Get BacMet database
mkdir ${PATH_TO_CONDA_ENVS}/MASTER-WP5g7/db/bacmet
wget -P http://bacmet.biomedicine.gu.se/download/BacMet2_EXP_database.fasta \
  ${PATH_TO_CONDA_ENVS}/MASTER-WP5g7/db/bacmet
wget -P wget http://bacmet.biomedicine.gu.se/download/BacMet2_EXP.753.mapping.txt \
  ${PATH_TO_CONDA_ENVS}/MASTER-WP5g7/db/bacmet
makeblastdb -in ${PATH_TO_CONDA_ENVS}/MASTER-WP5g7/db/bacmet/BacMet2_EXP_database.fasta \
  -parse_seqids -dbtype prot -hash_index -title BacMet2_EXP \
  -out ${PATH_TO_CONDA_ENVS}/MASTER-WP5g7/db/bacmet/BacMet2_EXP
```
> Modify the ${PATH_TO_CONDA_ENVS} and ${CPUS} variables accordingly.


<br>

## Running assembly-based pipeline

### Running from MAGs:

```
tormes --metadata metadata.txt --output ${OUTPUT_DIRECTORY} \
--threads ${CPUS} --no_mlst --no_pangenome --only_gene_prediction \
--prodigal_options “-p single” \
--gene_min_cov 80 --gene_min_id 80 \
--custom_genes_db "plasmidfinder"
```

### Running from whole assembly:

```
tormes --metadata metadata.txt --output ${OUTPUT_DIRECTORY} \
--threads ${CPUS} --no_mlst --no_pangenome --only_gene_prediction \
--prodigal_options “-p meta” \
--gene_min_cov 80 --gene_min_id 80 \
--custom_genes_db "plasmidfinder"
```

<br>

The parameters ```--gene_min_cov``` and ```--gene_min_id``` affect the screening of AMR and virulence genes and can be adjusted depending on the scope of the study.

<br>

The software to run here, [TORMES](https://github.com/nmquijada/tormes), requires a metadata file containing the name of the samples and the path to their fasta files. This metadata file has to be parsed with the ```-m/--metadata``` flag an it requires a specific organization:

- Columns must be tab separated.
- First column must be called ‘Samples’ and harbor samples’ names (avoid special characters and names composed only by numbers).
- Second column must be called ‘Read1’. As TORMES will be used with MAGs/contigs, This field has to contain the word "GENOME" (beware the capital letters!).
- Third column must be called ‘Read2’ and harbor the path to the genome/MAG/contig (in FASTA format)
- Fourth (and so on) columns are optional. The information included here is not needed for TORMES to work but will be included in the interactive report. You can add as many description columns as needed (information such as the sample where the MAG was extracted, MASTER partner, company, etc.).

Further information regarding the generation of this metadata file can be found [here](https://github.com/nmquijada/tormes#obligatory-options).

<br>

### Working with BacMet database (will be included in the next TORMES release)
```
# Use the *.faa files generated after gene prediction with Prodigal within TORMES
mkdir ${OUTPUT_DIRECTORY}/bacmet
blastp -query ${OUTPUT_DIRECTORY}/annotation/genome/genome.faa \
  -db ${PATH_TO_CONDA_ENVS}/MASTER-WP5g7/db/bacmet/BacMet2_EXP \
  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen' \
  -evalue 1e-25 -num_threads 60 -out ${OUTPUT_DIRECTORY}/bacmet/genome.blastout
# Calculate query coverage
cat ${OUTPUT_DIRECTORY}/bacmet/genome.blastout | awk -v OFS="\t" -F "\t" '{print $0, $14=($4-$6)*100/$13}' \
  > ${OUTPUT_DIRECTORY}/bacmet/genome.blast.bacmet.raw.txt
sed -i "1iqseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tslen\tqcov" \
  ${OUTPUT_DIRECTORY}/bacmet/genome.blast.bacmet.raw.txt

# Select those hits with more than 60% coverage and identity and just one hit per gene (based on e value)
tail -n+2 ${OUTPUT_DIRECTORY}/bacmet/genome.blast.bacmet.raw.txt | awk -v OFS="\t" -F "\t" '($3 > 60)' | \
  awk -v OFS="\t" -F "\t" '($14 > 60)' > ${OUTPUT_DIRECTORY}/bacmet/genome.blast.bacmet.raw.tmp
for i in $(cut -f 1 ${OUTPUT_DIRECTORY}/bacmet/genome.blast.bacmet.raw.tmp | sort -u); do
  grep "$i" ${OUTPUT_DIRECTORY}/bacmet/genome.blast.bacmet.raw.tmp | head -n 1 \ 
    >> ${OUTPUT_DIRECTORY}/bacmet/genome.blast.bacmet.min_cov_id-60.txt
done
sed -i "1i$(head -n 1 ${OUTPUT_DIRECTORY}/bacmet/genome.blast.bacmet.raw.txt)" \
  ${OUTPUT_DIRECTORY}/bacmet/genome.blast.bacmet.min_cov_id-60.txt

## Cleaning the house
rm -f ${OUTPUT_DIRECTORY}/bacmet/genome.blastout ${OUTPUT_DIRECTORY}/bacmet/genome.blast.bacmet.raw.tmp
```
> Modify "genome.faa" by either the MAGs or contigs files' names.

<br>

## Output from the assembly-based pipeline

### Annotation gff files for WP5g6 usage

One .gff, .faa and .nt file per MAG/contigs file will be created and can be found in:

```${OUTPUT_DIRECTORY}/gene_prediction/```

The gff files will be used for the metagenome annotation following WP5 Group 6 pipeline.

### AMR and virulence gene screening

The antimicrobial resistance genes screening results (against the ResFinder database) for each MAG/contigs will be found in: 

```${OUTPUT_DIRECTORY}/antibiotic_resistance_genes/resfinder/```

The antimicrobial resistance genes screening results (against the CARD database) for each MAG/contig will be found in:

```${OUTPUT_DIRECTORY}/antibiotic_resistance_genes/card/```

The virulence genes screening results (against the VFDB) for each MAG/contigs will be found in:

```${OUTPUT_DIRECTORY}/virulence_genes/```

The screening against BacMet (antibacterial biocide and metal resistance genes database) will be found in:

`${OUTPUT_DIRECTORY}/bacmet/`

<br>

## Running assembly-free pipeline

```
bowtie2 -x ${PATH_TO_CONDA_ENVS}/MASTER-WP5g7/db/resfinder/resfinder.bowtie \
-1 ${QUALITY_FILTERED_READS}_R1.fastq.gz -2 ${QUALITY_FILTERED_READS}_R2.fastq.gz \
-S ${OUTPUT_DIRECTORY}/${SAMPLE}.sam -p ${CPUS} --end-to-end \
--very-sensitive 2>>${OUTPUT_DIRECTORY}/${SAMPLE}-mapping-stats.txt

bowtie2 -x ${PATH_TO_CONDA_ENVS}/MASTER-WP5g7/db/card/card.bowtie \
-1 ${QUALITY_FILTERED_READS}_R1.fastq.gz -2 ${QUALITY_FILTERED_READS}_R2.fastq.gz \
-S ${OUTPUT_DIRECTORY}/${SAMPLE}.sam -p ${CPUS} --end-to-end \
--very-sensitive 2>>${OUTPUT_DIRECTORY}/${SAMPLE}-mapping-stats.txt

bowtie2 -x ${PATH_TO_CONDA_ENVS}/MASTER-WP5g7/db/vfdb/vfdb.bowtie \
-1 ${QUALITY_FILTERED_READS}_R1.fastq.gz -2 ${QUALITY_FILTERED_READS}_R2.fastq.gz \
-S ${OUTPUT_DIRECTORY}/${SAMPLE}.sam -p ${CPUS} --end-to-end \
--very-sensitive 2>>${OUTPUT_DIRECTORY}/${SAMPLE}-mapping-stats.txt
```

## Output from the assembly-free pipeline

A mapping file (.sam) will be created from each sample and will be found in the output directory (```${OUTPUT_DIRECTORY}/${SAMPLE}.sam```).
