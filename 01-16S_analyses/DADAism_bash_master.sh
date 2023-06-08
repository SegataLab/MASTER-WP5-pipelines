#!/bin/bash

# Author: Livio Antonielli (livio.antonielli@gmail.com) and Narciso Martin Quijada (nmartinquijada@gmail.com)

# Main software needed:
<<software
FASTQC
Bowtie2
Cutadapt
Metaxa2
MAFFT
BLAST
HMMER3
ITSx
fastx_len_sel.py script
R and Bioconductor
DADA2 and many more R packages: see R scripts for details about packages.
software

# Installation and usage
<<install
1. Install the environment first:
	conda env ceate -n MASTER-WP5g1 --file environment.yml

2. Activate the environment:
	conda activate MASTER-WP5g1

3. Copy the master bash script, the R scripts and the Python script in one directory and set the ${DADA_PATH}, accordingly.

4. Define the number of threads ${CPUS}.

5. Create a project directory and copy your Illumina paired-end FASTQ files in there (fastq.gz, fastq.tar.gz, or fastq.bz2 are also accepted).
In case your amplicons were generated from multiple forward and reverse primers, provide in the main directory a "fwd_primer.fasta" and a "rev_primer.fasta" file, respectively.

6. Launch the bash script inside the FASTQ directory.
install

# Set PATH (defined locally):
DADA_PATH=''

# ALIAS (defined locally, assuming scripts in DADA_PATH):
alias fastxLenSel="python ${DADA_PATH}/fastx_len_sel.py"
alias dadaPt1="Rscript ${DADA_PATH}/DADAism_pt.1.R"
alias dadaPt2="Rscript ${DADA_PATH}/DADAism_pt.2.R"
alias dadaPt3Bact="Rscript ${DADA_PATH}/DADAism_pt.3_bact.R"
alias dadaPt3Fung="Rscript ${DADA_PATH}/DADAism_pt.3_fung.R"

# Define parameters: 
CPUS=''


# Uncompressing files
for file in *.fastq*; do
	if [[ $file == *.fastq.bz2 ]]; then
		echo -e "\n ${file} fastq.bz2 uncompressing... \n"
		notify-send "Data pre-processing" "fastq bzp2 uncompressing and renaming" -t 10000
		bzip2 -cdk < $file > ${file//.bz2/}
	elif [[ $file == *.fastq.gz ]]; then
		echo -e "\n ${file} fastq.gz uncompressing... \n"
		notify-send "Data pre-processing" "fastq gz uncompressing and renaming" -t 10000
		gzip -cdk < $file > ${file//.gz/}
	elif [[ $file == *.fastq.tar.gz ]]; then
		echo -e "\n  ${file} fastq.tar.gz uncompressing... \n"
		notify-send "Data pre-processing" "fastq tar.gz uncompressing and renaming" -t 10000
		tar -czxvf < $file > ${file//.tar.gz/}
	elif [[ $file == *.fastq ]]; then
		echo -e "\n ${file} reads are in FASTQ format \n"
	else echo -e "\n ${file} format non supported, this will most likely not end up well... \n"
	fi
done


# Renaming samples keeping the first part of the name (assuming "_" as separator) 
ls *.fastq > fastq.list

while read R1; do
	read R2
	new_name=$(ls $R1 | awk -F "_" '{print $1}')
	mv $R1 ${new_name}_R1.fastq; mv $R2 ${new_name}_R2.fastq 
done < fastq.list


# Generating target directories
mkdir 1.reads 2.no_phix 3.stripped 4.filtered 5.output stats target_seq taxonomy_db 

# Raw read number stats: FASTQC
echo -e "\n==================== Counting reads... ====================\n"
mv *.fastq 1.reads/

for fq in 1.reads/*_R1.fastq; do                                               
	echo "$fq : `grep -c "^+$" $fq`"  
done > stats/1.raw_reads.counts

cat 1.reads/*_R1.fastq > 1.reads/raw_reads_R1.fastq
cat 1.reads/*_R2.fastq > 1.reads/raw_reads_R2.fastq
echo "Total read count : `cat 1.reads/raw_reads_R1.fastq | grep -c "^+$"`" >> stats/1.raw_reads.counts 
fastqc 1.reads/raw_reads_R1.fastq 1.reads/raw_reads_R2.fastq -o stats/
rm -rf 1.reads/raw_reads_R1.fastq 1.reads/raw_reads_R2.fastq


# PhiX filtering: Bowtie2
echo -e "\n==================== PhiX filtering... ====================\n"
ls 1.reads/*.fastq > 1.reads/fastq_raw.list

# Downloading PhiX sequence and building db
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/enterobacteria_phage_phix174_sensu_lato_uid14015/NC_001422.fna
bowtie2-build NC_001422.fna phix 

while read R1; do
	read R2
	echo $R1 $R2
	name=$(ls $R1 | awk -F "[/_]" '{print $2}')
	bowtie2 -x phix -1 $R1 -2 $R2 -t -p ${CPUS} --un-conc 2.no_phix/${name}.fastq -S 2.no_phix/${name}_contaminated_align.sam
	mv 2.no_phix/*.1.fastq 2.no_phix/${name}_R1.fastq
	mv 2.no_phix/*.2.fastq 2.no_phix/${name}_R2.fastq
done < 1.reads/fastq_raw.list

rm -rf 2.no_phix/*.sam *.fna *.bt2


# No-Phix read number stats: FASTQC
echo -e "\n==================== Counting phiX filtered reads... ====================\n"
for fq in 2.no_phix/*_R1.fastq; do                                               
	echo "$fq : `grep -c "^+$" $fq`"  
done > stats/2.nophix_reads.counts

cat 2.no_phix/*_R1.fastq > 2.no_phix/nophix_reads_R1.fastq
cat 2.no_phix/*_R2.fastq > 2.no_phix/nophix_reads_R2.fastq
echo "Total read count: `grep -c "^+$" 2.no_phix/nophix_reads_R1.fastq`" >> stats/2.nophix_reads.counts
fastqc  2.no_phix/nophix_reads_R1.fastq 2.no_phix/nophix_reads_R2.fastq -o stats/
rm -rf 2.no_phix/nophix_reads_R1.fastq 2.no_phix/nophix_reads_R2.fastq


# Primer stripping: Cutadapt
echo -e "\n==================== Primer stripping... ====================\n"
ls 2.no_phix/*.fastq > 2.no_phix/no_phix.list

read -p "Type '1' if you want to provide the primers manually (only one forward primer and one reverse primer) or type '2' if you have a list of primers: " trim_trick

if test "$trim_trick" = "1"
then
	#Cutadapt in manual mode 
	read -p "Enter your fwd primer sequence (5'-3'): " fwd
	read -p "Enter your rev primer sequence (5'-3'): " rev

	while read R1; do
		read R2
		cutadapt -j ${CPUS} -g $fwd -G $rev -e .2 --pair-filter=any --minimum-length 100 --discard-untrimmed -o 3.stripped/${R1##*/} -p 3.stripped/${R2##*/} $R1 $R2 
	done < 2.no_phix/no_phix.list
else
	#Cutadapt multi-primer mode
	while read R1; do
		read R2
		cutadapt -j ${CPUS} -g file:fwd_primer.fasta -G file:rev_primer.fasta -e .2 --pair-filter=any --minimum-length 100 --discard-untrimmed -o 3.stripped/${R1##*/} -p 3.stripped/${R1##*/} $R1 $R2
	done < 2.no_phix/no_phix.list 
fi


# Trimmed read number stats: FASTQC
echo -e "\n==================== Counting trimmed reads... ====================\n"
for fq in 3.stripped/*_R1.fastq; do                                               
	echo "$fq : `grep -c "^+$" $fq`"  
done > stats/3.stripped_reads.counts

cat 3.stripped/*_R1.fastq > 3.stripped/stripped_reads_R1.fastq
cat 3.stripped/*_R2.fastq > 3.stripped/stripped_reads_R2.fastq
echo "Total read count: `grep -c "^+$" 3.stripped/stripped_reads_R1.fastq`" >> stats/3.stripped_reads.counts
fastqc 3.stripped/stripped_reads_R1.fastq 3.stripped/stripped_reads_R2.fastq -o stats/
rm -rf 3.stripped/stripped_reads_R1.fastq 3.stripped/stripped_reads_R2.fastq


# Setting directory paths for R script
currentpath=$(pwd)
seqpath=$(readlink -f 3.stripped)


# Launching DADA2 script pt.1
echo -e "\n==================== R script pt.1... ====================\n"
dadaPt1 $currentpath $seqpath


# Moving diagnostic plots in stats folder
echo -e "\n==================== Please check diagnostic plots in stats folder... ====================\n"
mkdir stats/diagnostic_plots
mv *.png stats/*fastqc.html stats/*fastqc.zip stats/diagnostic_plots 

read -p "After diagnostic plot inspection, truncate R1 reads at this length: " truncLenR1
read -p "After diagnostic plot inspection, truncate R2 reads at this length: " truncLenR2


# Launching DADA2 script pt.2
echo -e "\n==================== R script pt.2... ====================\n"
filtpath=$(readlink -f 4.filtered)
dadaPt2 $currentpath $seqpath $truncLenR1 $truncLenR2 $filtpath

mv *.pdf stats/diagnostic_plots
mv *.png stats/diagnostic_plots
mv *.counts stats/4.filtered_reads.counts
mv *.fasta 5.output
mv seqtab*.txt 5.output
seqtab=$(readlink -f 5.output/seqtab_head_names.txt)


# Choosing taxonomic db (https://benjjneb.github.io/dada2/assign.html#species-assignment)
echo -e "\n==================== Downloading DB for taxonomic classification... ====================\n"
read -p "Choose your amplicon target between '16S rRNA' or 'fungal ITS' (type either 16S or ITS): " seq_type

if test "$seq_type" = "16S"
then
	# Download SILVA db v.138.1 and launch Metaxa2
	wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz
	wget https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz
	mv *.fa.gz taxonomy_db
	silva=$(readlink -f taxonomy_db/silva_nr99_v138.1_train_set.fa.gz)
	silvaSpec=$(readlink -f taxonomy_db/silva_species_assignment_v138.1.fa.gz)
	echo -e "\n==================== Metaxa2: 16s rRNA targeted extraction and verification... ====================\n"
	metaxa2 -i 5.output/seqs.fasta -o target_seq/target -t all --complement T --truncate T --plus T --cpu ${CPUS}
	sed -i '/>/ s/|.*//g' target_seq/target.extraction.fasta
	fastxLenSel target_seq/target.extraction.fasta 
	mv sel_seqs.fasta 5.output/seqs_targ.fasta 
else
	# Download UNITE db v.8.2 and launch ITSx
		cd taxonomy_db	
		wget https://files.plutof.ut.ee/public/orig/1D/B9/1DB95C8AC0A80108BECAF1162D761A8D379AF43E2A4295A3EF353DD1632B645B.gz
		mv *.gz sh_general_release_s_all_04.02.2020.fasta.gz
		gunzip sh_general_release_s_all_04.02.2020.fasta.gz
		sed -i 's/.*\(>.*\)/\1/' sh_general_release_s_all_04.02.2020.fasta
		cd ..
		unite=$(readlink -f taxonomy_db/*.fasta)
		read -p "Choose your amplicon target between ITS1 or ITS2: " seq_type_its

		if test "$seq_type_its" = "ITS1"
		then
			echo -e "\n==================== ITSx: ITS1 targeted extraction and verification... ====================\n"	 
			itsx -i 5.output/seqs.fasta -o target_seq/target --preserve F --only_full T --complement T -E 0.001 --plus T --cpu ${CPUS}
			sed -i '/>/ s/|.*//g' target_seq/target.ITS1.fasta
			fastxLenSel target_seq/target.ITS1.fasta
			mv sel_seqs.fasta 5.output/seqs_targ.fasta
		else
			echo -e "\n==================== ITSx: ITS2 targeted extraction and verification... ====================\n"
			itsx -i 5.output/seqs.fasta -o target_seq/target --preserve F --only_full T --complement T -E 0.001 --plus T --cpu ${CPUS}
			sed -i '/>/ s/|.*//g' target_seq/target.ITS2.fasta
			fastxLenSel target_seq/target.ITS2.fasta
			mv sel_seqs.fasta 5.output/seqs_targ.fasta
	fi
fi 


# Launching DADA2 script pt.3
echo -e "\n==================== R script pt.3... ====================\n"
asvSeqsTarg=$(readlink -f 5.output/seqs_targ.fasta)

if test "$seq_type" = "16S"
then
	# Launching DADA2 script pt.3 for bacteria
	dadaPt3Bact $currentpath $seqtab $asvSeqsTarg $silva $silvaSpec
else
	if test "$seq_type" = "ITS"
	then
		# Launching DADA2 script pt.3 for fungi
		dadaPt3Fung $currentpath $seqtab $asvSeqsTarg $unite	 
	fi
fi 


# Moving goodies in stats and output folders
mv *.counts stats
mv *.txt 5.output


# Creating 2 FASTA files out of the taxon_seq_table, with a simple header and also the taxonomy, respectively 
while read line; do 
	echo $line | cut -d" " -f1,9 | sed 's\"\\g' | awk -F" " '{ print ">"$1"\n"$2 }' >> 5.output/asv_seqs.fasta
done < <(sed '1d' 5.output/taxon_seq_table.txt)

while read line; do 
	echo $line | sed 's\"\\g' | awk -F" " '{ print ">"$1"\t"$2";"$3";"$4";"$5";"$6";"$7";"$8";""\n"$9 }' >> 5.output/asv_seqs_tax.fasta
done < <(sed '1d' 5.output/taxon_seq_table.txt)


# Data processing End 
echo -e "\n==================== Data are ready for decontam and further statistics... ====================\n"
