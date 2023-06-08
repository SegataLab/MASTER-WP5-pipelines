#!/bin/bash

### The code assume that inside a master folder with absolute path pathReads=/path/${dataset_name}/reads there is a folder for each sample (and named after the sample), inside which there are the files with the reads. The files are in fastq format and zipped with respective name ${samplename}_R1.fastq.bz2, ${samplename}_R2.fastq.bz2, ${samplename}_UN.fastq.bz2
#i.e. /path/${dataset_name}/reads/${samplename}/${samplename}_R1.fastq.bz2


samples=`ls -1 /file | tr '\n' ' ' `

## STEP 1: assembly of reads in contigs using megahit
NOJ= #NumberOfJobs: insert reasonable number of parallel jobs according to hardware available
parallel -j${NOJ} 'run_single_assembly.sh ${dataset_name} {}' ::: ${samples}

#STEP 2: filter contigs according to length
check python version and packages needed
NOJ= #NumberOfJobs: insert reasonable number of parallel jobs according to hardware available
parallel -j${NOJ} 'place="${dataset_name}/assemblies_output"; python filter_contigs.py ${place}/{}/final.contigs.fa -n {} -l ${minimum_read_length} > ${place}/{}/filtered_contigs.fa' ::: ${samples}

#STEP 3: align filtered contigs agianst original reads in order ...
NOJ= #NumberOfJobs: insert reasonable number of parallel jobs according to hardware available
NOTh= #NumberOfThreads: insert reasonable number of parallel jobs according to hardware available
parallel -j${NOJ} 'place="${dataset_name}/assemblies_output"; pathReads="/path/${dataset_name}/reads"; bowtie2-build ${place}/{}/filtered_contigs.fa ${place}/{}/{}; if [ -f "${pathReads}/{}/{}_R1.fastq.bz2"  ]; then bowtie2 -x ${place}/{}/{} -1 ${pathReads}/{}/{}_R1.fastq.bz2 -2 ${pathReads}/{}/{}_R2.fastq.bz2 -S - --very-sensitive-local --no-unal -p ${NOTh} | samtools view -bS - > ${place}/{}/{}.unsorted.bam; fi; if [ -f "${pathReads}/{}/{}_UN.fastq.bz2" ]; then bowtie2 -x ${place}/{}/{} -U ${pathReads}/{}/{}_UN.fastq.bz2 -S - --very-sensitive-local --no-unal -p ${NOTh} | samtools view -bS - >> ${place}/{}/{}.unsorted.bam; fi; samtools sort ${place}/{}/{}.unsorted.bam -o ${place}/{}/{}.bam -m 8G; samtools index ${place}/{}/{}.bam ${place}/{}/{}.bam.bai' ::: ${samples}

#STEP 4: use metabat to find contig depths
NOJ= #NumberOfJobs: insert reasonable number of parallel jobs according to hardware available
parallel -j${NOJ} 'place="${dataset_name}/assemblies_output"; jgi_summarize_bam_contig_depths --outputDepth ${place}/{}/depth.txt ${place}/{}/{}.bam' ::: ${samples}

#STEP 5: use metabat to compact contigs into bins/putative MAGs
NOJ= #NumberOfJobs: insert reasonable number of parallel jobs according to hardware available
NOTh= #NumberOfThreads: insert reasonable number of parallel jobs according to hardware available
parallel -j${NOJ} 'place="${dataset_name}/assemblies_output"; mkdir -p ${place}/{}/bins_dir; metabat2 -m 1500 -t ${NOTh} --unbinned --seed 0 -i ${place}/{}/filtered_contigs.fa -a ${place}/{}/depth.txt -o ${place}/{}/bins_dir/bin' ::: ${samples}

#STEP 6: use checkm to verify completeness and contamination of the reconstructed bin and choose which are adequate MAGs
NOJ= #NumberOfJobs: insert reasonable number of parallel jobs according to hardware available
NOTh= #NumberOfThreads: insert reasonable number of parallel jobs according to hardware available

export PATH=path/to/tools/hmmer-3.1b2/binaries:$PATH
export PATH=path/to/tools/prodigal-2.6.3:$PATH
export PATH=path/to/tools/pplacer-1.1.alpha19:$PATH

parallel -j5 'place="${dataset_name}/assemblies_output"; checkm lineage_wf -t ${NOTh} -x fa ${place}/{}/bins_dir ${place}/{}/bins_dir/checkm > ${place}/{}/bins_dir/{}.log' ::: ${samples}

#final (optional): copy all MAGs surpassing quality thresholds in a folder and create a summary file
./summarize_checkm.sh ${dataset_name}
