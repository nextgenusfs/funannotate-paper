#!/usr/bin/env python

#This script removes fasta files from list in file

import sys
from Bio.SeqIO.FastaIO import SimpleFastaParser

def softwrap(string, every=80):
    lines = []
    for i in xrange(0, len(string), every):
        lines.append(string[i:i+every])
    return '\n'.join(lines)

if len(sys.argv) < 3:
    print "Usage: fasta_remove.py input.fasta list2remove.txt"
    sys.exit(1)

#get list of names from file
remove = []
with open(sys.argv[2], 'rU') as filein:
	for line in filein:
		line = line.replace('\n', '').rstrip()
		if not line in remove:
			remove.append(line)

with open(sys.argv[1], 'rU') as input:
	for header, sequence in SimpleFastaParser(input):
		if not header in remove:
			sys.stdout.write('>{:}\n{:}\n'.format(header, softwrap(sequence)))



