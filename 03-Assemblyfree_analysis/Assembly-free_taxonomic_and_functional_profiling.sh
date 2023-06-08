

#This pipeline assumes that trimmed metagenomics reads in "fastq.bz2" format are placed in a folder namerd "fastq" to be submitted as inputs in strain tracking pipeline

#Create conda environment for Metaphlan3 and Humann3 software, included in biobakery3 library
#Ommit this step if biobakery3 is already installed
#conda create --name biobakery3 python=3.7
#conda activate biobakery3
#conda install humann -c biobakery3
#conda install bbmap -c biobakery3

#Create conda environment for kraken2 and bracken software
#Ommit this step if biobakery3 is already installed
#conda create --yes -n kraken kraken2 bracken
#Build full kraken database
#kraken2-build --standard --threads 24 -db $DATABASE_PATH

#Create conda environment for Superfocus software
#Ommit this step if Superfocus is already installed
#conda create --name Superfocus python=3.7
#conda activate Superfocus
#git clone https://github.com/metageni/SUPER-FOCUS.git
#cd SUPER-FOCUS && python setup.py install

#Source conda to activate conda environments
source /home/raulc/miniconda3/etc/profile.d/conda.sh

#Select the number of cores available in the workstation
n_cores=24

#METAPHLAN 3.0 ANALYSIS

#Activate conda environment for metaphlan 
conda activate biobakery3
#Create folders to store SAM, bowtie2 and profiling results
mkdir -p metaphlan3_sams
mkdir -p metaphlan3_bowtie2
mkdir -p metaphlan3_profiles

kraken2_db=$PATH
	
#Launch MetaPhlAn 3.0
for reads in fastq/*R1.fastq.bz2
do
    echo "Running MetaPhlAn on ${reads}"
    base_name=$(basename ${reads})
    fname="${base_name%R1.fastq.bz2}"
    metaphlan fastq/${fname}R1.fastq.bz2,fastq/${fname}R2.fastq.bz2 --input_type fastq -s metaphlan3_sams/${fname}sam.bz2 --bowtie2out metaphlan3_bowtie2/${fname}bowtie2.bz2 -o metaphlan3_profiles/${fname}profiled.tsv --nproc "${n_cores}"
    merge_metaphlan_tables.py metaphlan3_profiles/*.tsv -o metaphlan3_profiles/metaphlan.tsv
done

conda deactivate



#KRAKEN2 ANALYSIS

#Activate conda environment for kraken2 and bracken 
conda activate kraken

#Launch Kraken2/Bracken
for reads in fastq/*R1.fastq.bz2
do
    echo "Running KRAKEN2 on ${reads}"
    base_name=$(basename ${reads})
    fname="${base_name%R1.fastq.bz2}"
	kraken2 --memory-mapping --db "${kraken2_db}" --threads "${n_cores}" --report ${fname}taxonomy_report.txt --output ${fname}taxonomy.txt --bzip2-compressed --paired fastq/${fname}R1.fastq.bz2 fastq/${fname}R2.fastq.bz2
	bracken -d "${kraken2_db}" -i ${fname}taxonomy_report.txt -l S -r 100 -t "${n_cores}" -o ${fname}Species_bracken_output
	bracken -d "${kraken2_db}" -i ${fname}taxonomy_report.txt -l G -r 100 -t "${n_cores}" -o ${fname}Genera_bracken_output
	bracken -d "${kraken2_db}" -i ${fname}taxonomy_report.txt -l F -r 100 -t "${n_cores}" -o ${fname}Family_bracken_output
	bracken -d "${kraken2_db}" -i ${fname}taxonomy_report.txt -l P -r 100 -t "${n_cores}" -o ${fname}Phylum_bracken_output
	
done

#Store kraken2 and bracken2 results in a different folder
mkdir kraken2_bracken_results
mv *taxonomy* kraken2_bracken_results/
mv *bracken* kraken2_bracken_results/

combine_bracken_outputs.py --files kraken2_bracken_results/*_Phylum_bracken_output -o combined_phylum_mpa.csv
combine_bracken_outputs.py --files kraken2_bracken_results/*_Family_bracken_output -o combined_family_mpa.csv
combine_bracken_outputs.py --files kraken2_bracken_results/*_Genera_bracken_output -o combined_genus_mpa.csv
combine_bracken_outputs.py --files kraken2_bracken_results/*_Species_bracken_output -o combined_species_mpa.csv

conda deactivate



#HUMANN 3.0 ANALYSIS

#Activate conda environment for humann 
conda activate biobakery3

for reads in fastq/*R1.fastq.bz2
do
    echo "Running HUMANN on ${reads}"
    base_name=$(basename ${reads})
    fname="${base_name%R1.fastq.bz2}"
	reformat.sh in=fastq/${fname}R1.fastq.bz2 in2=fastq/${fname}R2.fastq.bz2 out=${fname}_interleaved.fastq.gz
	humann -i ${fname}interleaved.fastq.gz -o ${fname}humann3out --threads "${n_cores}"  

	mkdir -p Humann3outputs/Genefamilies
	mkdir -p Humann3outputs/Pathabundance
	mkdir -p Humann3outputs/Pathcoverage
	mv ${fname}_humann3out/${fname}_genefamilies.tsv Humann3outputs/Genefamilies/${fname}.tsv
	mv ${fname}_humann3out/${fname}_pathabundance.tsv Humann3outputs/Pathabundance/${fname}.tsv
	mv ${fname}_humann3out/${fname}_pathcoverage.tsv Humann3outputs/Pathcoverage/${fname}.tsv
	rm -r "$i"_humann3out/

done

#Join all gene family and pathway abudance files
humann_join_tables -i Humann3outputs/Pathabundance/ -o Humann3outputs/pathabundance.tsv -s --file_name '.tsv'
humann_join_tables -i Humann3outputs/Pathcoverage/ -o Humann3outputs/pathcoverage.tsv -s --file_name '.tsv'
humann_join_tables -i Humann3outputs/Genefamilies/ -o Humann3outputs/genefamilies.tsv -s --file_name '.tsv'

conda deactivate



#SUPERFOCUS ANALYSIS

#Activate conda environment for humann 
conda activate Superfocus

#Create imput folders and descompress
mkdir -p FastqSup
cp fastq/*.bz2 FastqSup/
bzip2 -d FastqSup/*.bz2

#Launch SUPERFOCUS
superfocus -q FastqSup -a diamond -t "${n_cores}" -dir superfocus_result

rm -R FastqSup/

conda deactivate



