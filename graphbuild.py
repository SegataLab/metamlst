#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import sqlite3
import time,os.path
import argparse

conn = sqlite3.connect('bsb.db')
conn.row_factory = sqlite3.Row
c = conn.cursor()
parser = argparse.ArgumentParser()
parser.add_argument("organism", help="bacterium name (e.g. saureus)")
parser.add_argument("infile", help="Queue these sequences")
parser.add_argument("-a","--add", help="Use this unique file to queue")
parser.add_argument("-e","--exclude", help="Exclude these genes (comma separated)")
args=parser.parse_args()

if args.add: outName = args.add
else: outName =  args.organism+'_'+str(int(time.time()))+'.fasta'

if (args.add and not os.path.exists(args.add)) or not args.add:
	orderedGenes = sorted([row['geneName'] for row in c.execute("SELECT geneName FROM genes WHERE bacterium = ?",(args.organism,))], key = lambda x: x[0])
	
	profiles = {}
	if args.exclude: excludeList = "','".join([str(a) for a in (args.exclude).split(',')])
	else: excludeList = ''

	for row in c.execute("SELECT * FROM profiles,alleles WHERE alleleCode = alleles.recID AND gene NOT IN ('"+excludeList+"') AND profiles.bacterium = ?", (args.organism,)):
		if row['profileCode'] not in profiles: profiles[row['profileCode']] = {}
		profiles[row['profileCode']][row['gene']] = row['alignedSequence']
	print "SELECT * FROM profiles,alleles WHERE alleleCode = alleles.recID AND gene NOT IN ('"+excludeList+"') AND profiles.bacterium = ?", (args.organism,)
	recs = []
	for profile,genes in profiles.items():
		sequ = ''
		for gene,sequence in sorted(genes.items(), key = lambda x: x[0]):
			sequ = sequ+sequence
		
		recs.append(SeqRecord(Seq(sequ, IUPAC.unambiguous_dna), id = args.organism+'_'+str(profile), description = ''))
		
	SeqIO.write(recs, outName, "fasta")

dfil = open(outName,'a')

for record in SeqIO.parse(args.infile,'fasta'):
	dfil.write('>'+record.id+'\n')
	dfil.write(str(record.seq)+'\n')
dfil.close()
