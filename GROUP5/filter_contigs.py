#!/usr/bin/env python

import glob
import os
import sys
import argparse as ap
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def read_params(args):
	parser = ap.ArgumentParser(description='Filter contigs')
	arg = parser.add_argument
	arg( 'inp_f', metavar='INPUT_FILE', nargs='?', default=sys.stdin, type=str, help="the input fasta file [stdin if not present]")	
	arg( 'out_f', metavar='OUTPUT_FILE', nargs='?', default=None, type=str, help="the output file [stdout if not present]")	
	arg( '-l','--minimum_read_length', default=1000, type=int, help="minimum read length\n")
	arg( '-n','--name_to_append', default=None, type=str)
	return vars(parser.parse_args())


if __name__ == '__main__':
	par = read_params(sys.argv)

	if par['out_f']:
		fid = open(par['out_f'],'w')
	else:
		fid = sys.stdout


        #c = 0
        ##for record in SeqIO.parse(par['inp_f'], 'fasta'): #print len(record.seq), 
        #    if len(record.seq) > par['minimum_read_length']: c += 1
        #exit()

        #print c, ' is count'
        #print len([SeqRecord(record.seq, '_'.join([par['name_to_append'], record.id, record.description]), description='')
        #             for record in SeqIO.parse(par['inp_f'], 'fasta')
        #                 if len(record.seq) > par['minimum_read_length']]), ' is cut'

	SeqIO.write([SeqRecord(record.seq, '_'.join([par['name_to_append']\
                             , record.description]), description='') ## record.id,
                     for record in SeqIO.parse(par['inp_f'], 'fasta')
                     if len(record.seq) > par['minimum_read_length']],
                    fid, 'fasta')

	if par['out_f']:
		fid.close()
