#STRAIN TRACKING AND PROFILING METHOD

#This pipeline assumes that the names of microbial species of interest are in a txt file named "species_of_interest.txt" (each line corresponds to the name of one species, i. e. Eubacterium_rectale, Bacteroides_caccae)

#This pipeline assumes that trimmed metagenomics reads in "fastq.bz2" format are placed in a folder namerd "fastq" to be submitted as inputs in strain tracking pipeline

#This pipeline assumes that reference genomes (named as "ref_Genus_species") in "fna" format are placed in an additional folder named "reference_genomes" 
#Bacteroides_caccae reference genome: https://www.ebi.ac.uk/ena/browser/view/GCA_016726305?show=chromosomes
#Eubacterium_rectale reference genome: https://www.ebi.ac.uk/ena/browser/view/CP001107.1

#Create conda environment for strainphlan software, included in biobakery3 library
#Ommit this step if biobakery3 is already installed
#conda create --name biobakery3 python=3.7
#conda activate biobakery3
#conda install humann -c biobakery

#Create a new conda environment for GraPhlAn software to avoid python incompatibility issues
#Ommit this step if GraPhlAn is already installed
#conda create --name graphlan_env python=2.7
#Activate new conda environment for GraPhlAn software and installl GraPhlAn and dendropy modules
#Ommit this step if GraPhlAn is already installed
#conda activate graphlan_env
#pip install graphlan
#pip install dendropy
#If errors persist: place "plot_tree_graphlan.py" from the official repository to "bin" folder in GraPhlAn environment

#Create conda environment for panphlan software
#Ommit this step if biobakery3 is already installed
#conda create --name panphlan_env python=3.7
#conda activate panphlan_env
#conda install -c bioconda panphlan

#Create conda environment for kraken2 and bracken software
#Ommit this step if biobakery3 is already installed
#conda create --yes -n kraken kraken2 bracken
#Build full kraken database
#kraken2-build --standard --threads 24 -db $DATABASE_PATH

#Source conda to activate conda environments
source /home/usuario/anaconda3/etc/profile.d/conda.sh

#Select the number of cores available in the workstation
n_cores=24

#Select kraken2 database path
kraken2_db=/home/usuario/anaconda3/envs/tormes_env/db/minikraken-DB/minikraken_8GB_20200312

#Iterate the text file containing the names of microbial species of interest to perform strain tracking analyses
while read -r species; do (
	
	#STRAINPHLAN 3.0 ANALYSIS

	#Activate conda environment for strainphlan 
	conda activate biobakery3
	#Create folders to store SAM files, consensus and database markers as well as bowtie2 and profiling results
	mkdir -p "${species}"_sams
	mkdir -p "${species}"_bowtie2
	mkdir -p "${species}"_profiles
	mkdir -p "${species}"_consensus_markers
	mkdir -p "${species}"_db_markers
	mkdir -p "${species}"_strainphlan_output
	
	#Generate SAM files for MASTER samples by mapping them againts MetaPhlAn 3.0 marker database
	for reads in fastq/*fastq*
	do
	    echo "Running MetaPhlAn on ${reads}"
	    base_name=$(basename ${reads})
	    metaphlan ${reads} --input_type fastq -s "${species}"_sams/"${base_name}".sam.bz2 --bowtie2out "${species}"_bowtie2/"${base_name}".bowtie2.bz2 -o "${species}"_profiles/"${base_name}"_profiled.tsv
	done
	
	#Generate consensus-marker files from SAM files: reconstruct and store in a pickle file (pkl) all species strains present in samples using the StrainPhlAn 3.0 script “sample2markers.py”
	sample2markers.py -i "${species}"_sams/*.sam.bz2 -o "${species}"_consensus_markers -n "${n_cores}"	
	
	#Extract specific markers of microbial species of interest from ChocoPhlAn database using the StrainPhlAn 3.0 script “extract_markers.py”
	#This example shows how to exctract marker genes from Bacteroides caccae using ChocoPhlAn database:
	extract_markers.py -c s__"${species}" -o "${species}"_db_markers/
	
	#Clean seq-id names from reference genomes
	sed -i -e 's/|/_/g' reference_genomes/ref_"${species}".fasta
	
	#Filter selected clade markers based on their presence in the consensus marker files (where all strains reconstructed are stored). As a result, filtered markers are be obtained and 		StrainPhlAn 3.0 in combination with PhyloPhlAn 3.0 software are  used to produce a multiple sequence alignment (MSA) and phylogenetic trees
	strainphlan -s "${species}"_consensus_markers/*.pkl -m "${species}"_db_markers/s__"${species}".fna -r reference_genomes/ref_"${species}".fasta -o "${species}"_strainphlan_output -n "${n_cores}" -c s__"${species}" --phylophlan_mode accurate --mutation_rates
	
	#Metadata files can be added to phylogenetic trees using the StrainPhlAn 3.0 script "add_metadata_tree.py"
	add_metadata_tree.py -t "${species}"_strainphlan_output/RAxML_bestTree.s__"${species}".StrainPhlAn3.tre -f fastq/metadata.txt -m subjectID --string_to_remove .fastq.bz2

	#Activate conda environment for GraPhlAn software to plot the phylogenetic tree
	conda activate graphlan_env

	#Plot phylogenetic trees using the StrainPhlAn 3.0 script "plot_tree_graphlan.py"
	plot_tree_graphlan.py -t "${species}"_strainphlan_output/RAxML_bestTree.s__"${species}".StrainPhlAn3.tre.metadata -m subjectID

	#Activate strainphlan conda environment for strain transmission analysis
	conda activate biobakery3
	
	#Copy tree file to current folder in order to perform strain transmission analysis
	cp "${species}"_strainphlan_output/RAxML_bestTree.s__"${species}".StrainPhlAn3.tre $PWD
	
	#Execute StrainPhlAn 3.0 script “strain-transmission.py” to perform strain-tracking analysis
	strain_transmission.py -t RAxML_bestTree.s__"${species}".StrainPhlAn3.tre -m fastq/metadata.tsv --output_dir "${species}"_strainphlan_output --threshold 0.1
	
	#Delete tree file copied to perform strain transmission analysis
	rm -rf *.tre	
	
	#Metadata file should contain information about subjects, timepoints and samples
	#relation = line[2]
	#subject = line[1]
	#timepoint = line[3]
	#sample = line[0]




	#PANHPHLAN 3.0 ANALYSIS
	
	#Activate conda panphlan environment 
	conda activate panphlan_env
	
	#Download reference pangenomes
	panphlan_download_pangenome.py -i "${species}" -o "${species}"
	
	#Clean reference pangenomes to avoid formatting errors
	panphlan_clean_pangenome.py --species "${species}" --pangenome "${species}"/"${species}"
	
	#Map samples against reference pangenome
	for reads in fastq/*fastq*; do (
	    base_name=$(basename $reads)
	    cd "${species}"/"${species}"
	    panphlan_map.py -i ../../samples_fastq/${base_name} --indexes "${species}" -p "${species}"_pangenome.tsv -o ../../"${species}"_map_results/${basename%.*}_"${species}".tsv
	    cd ../../
	); done

	#Strain profiling analysis
	panphlan_profiling.py -i "${species}"_map_results/ --o_matrix result_profile_"${species}".tsv -p "${species}"/"${species}"/"${species}"_pangenome.tsv --add_ref 
	
	
	
	
	#KRAKEN2 ANALYSIS
	
	#Activate conda environment for kraken2 and bracken 
	conda activate kraken

	#Iterate samples reads
	for reads in fastq/*fastq*; do (
		base_name=$(basename $reads)
		
		#Kraken2 taxonomic classification
		kraken2 --memory-mapping --db "${kraken2_db}" --threads "${n_cores}" --output "${basename%.*}"_taxonomy.txt --report "${basename%.*}"_taxonomy_report.txt $reads
		
		#Post-process kraken2 results using bracken at strain level (S1)
		bracken -d "${kraken2_db}" -i "${basename%.*}"_taxonomy_report.txt -l S1 -o "${basename%.*}"_bracken_output
		
	); done
	
	#Store kraken2 and bracken2 results in a different folder
	mkdir kraken2_bracken_results
	mv *taxonomy* kraken2_bracken_results/
	mv *bracken* kraken2_bracken_results/
	
); done < species_of_interest.txt

