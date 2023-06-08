#!/usr/bin/python3

#Livio Antonielli (livio.antonielli@gmail.com)

"""
The script calculates basic statistics on read length and number of reads when a FASTA or FASTQ file (also FASTQ gunzip compressed) is provided. 
Then, a minimum and maximum read length number must be provided to filter the reads accordingly.
Use as script.py seq.fasta 
"""

import gzip
import sys
from Bio import SeqIO
import statistics

seq_file = sys.argv[1]
records = 0

if seq_file.endswith('.fastq.gz'):
    records = list(SeqIO.parse(gzip.open(seq_file, 'rt'), "fastq"))
elif seq_file.endswith('.fastq'):
    records = list(SeqIO.parse(str(seq_file), "fastq"))
elif seq_file.endswith('.fa'):
    records = list(SeqIO.parse(str(seq_file), "fasta"))    
elif seq_file.endswith('.fna'):
    records = list(SeqIO.parse(str(seq_file), "fasta"))     
elif seq_file.endswith('.fasta'):
    records = list(SeqIO.parse(str(seq_file), "fasta"))       
else:
    print("Sequence file extension not supported!")
    quit()

seq_nums = len(records)
print("Total reads: %i" % seq_nums)

seq_lengths = [len(record) for record in records]
seq_lengths_mean = statistics.mean(seq_lengths)
seq_lengths_median = statistics.median(seq_lengths)
seq_lengths_min = min(seq_lengths)
seq_lengths_max = max(seq_lengths)

seq_lengths_mean_records = [] # Setup an empty list
for record in records:
    if len(record.seq) == int(seq_lengths_mean):
        seq_lengths_mean_records.append(record)
        
seq_lengths_median_records = [] 
for record in records:
    if len(record.seq) == int(seq_lengths_median):
        seq_lengths_median_records.append(record)        

seq_lengths_min_records = [] 
for record in records:
    if len(record.seq) == int(seq_lengths_min):
        seq_lengths_min_records.append(record)
        
seq_lengths_max_records = []
for record in records:
    if len(record.seq) == int(seq_lengths_max):
        seq_lengths_max_records.append(record)        

print("Mean: %.2f bp (%i reads, %.2f%%)" % (seq_lengths_mean, len(seq_lengths_mean_records), len(seq_lengths_mean_records) / seq_nums * 100))
print("Median: %.2f bp (%i, %.2f%%)" % (seq_lengths_median, len(seq_lengths_median_records), len(seq_lengths_median_records) / seq_nums * 100))
print("Min: %i bp (%i, %.2f%%)" % (seq_lengths_min, len(seq_lengths_min_records), len(seq_lengths_min_records) / seq_nums * 100))
print("Max: %i bp (%i, %.2f%%)" % (seq_lengths_max, len(seq_lengths_max_records), len(seq_lengths_max_records) / seq_nums * 100))
print()

for seq_length in sorted(set(seq_lengths)):
    match_length_records = []
    for record in records:
        if len(record.seq) == seq_length:
            match_length_records.append(record)
    print("%i bp\t%i reads\t(%.2f%%)" % (seq_length, len(match_length_records), len(match_length_records) / seq_nums * 100))
print()
            
seq_len_min = input('Enter minimum acceptable read length: ')
seq_len_max = input('Enter maximum acceptable read length: ')
print()

range_lengths_records = []
for record in records:
    if len(record.seq) in range(int(seq_len_min), int(seq_len_max)+1):
        range_lengths_records.append(record)

print("Found %i reads with length between %i and %i bp." % (len(range_lengths_records), int(seq_len_min), int(seq_len_max)))

SeqIO.write(range_lengths_records, "sel_seqs.fasta", "fasta")
