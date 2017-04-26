#!/usr/bin/env python

try:
	import argparse
	import math
	import sqlite3
	import subprocess
	import time
	import gc
	import os
	import sys
	from metaMLST_functions import *
except ImportError as e:
	print "Error while importing python modules! Remember that this script requires: sys,os,subprocess,sqlite3,argparse,re"
	sys.exit(1)

try:
	from Bio import SeqIO
	from Bio.Seq import Seq
	from Bio.SeqRecord import SeqRecord
	from Bio.Alphabet import IUPAC
except ImportError as e:
	metamlst_print("Failed in importing Biopython. Please check Biopython is installed properly on your system!",'FAIL',bcolors.FAIL)
	sys.exit(1)
 
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		description='Reconstruct the MLST loci from a BAMFILE aligned to the reference MLST loci')

parser.add_argument("BAMFILE", help="BowTie2 BAM file containing the alignments",nargs='?')
parser.add_argument("-o", metavar="OUTPUT FOLDER", help="Output Folder (default: ./out)", default='./out')
parser.add_argument("-d", metavar="DB PATH", help="MetaMLST SQLite Database File (created with metaMLST-index)", default=os.path.abspath(os.path.dirname(__file__))+'/metamlstDB_2017.db')

parser.add_argument("--filter", metavar="species1,species2...", help="Filter for specific set of organisms only (METAMLST-KEYs, comma separated. Use metaMLST-index.py --listspecies to get MLST keys)")
parser.add_argument("--penalty", metavar="PENALTY", help="MetaMLST penaty for under-represented alleles", default=100, type=int)
parser.add_argument("--minscore", metavar="MINSCORE", help="Minimum alignment score for each alignment to be considered valid", default=80, type=int)
parser.add_argument("--max_xM", metavar="XM", help="Maximum SNPs rate for each alignment to be considered valid (BowTie2s XM value)", default=5, type=int)
parser.add_argument("--min_read_len", metavar="LENGTH", help="Minimum BowTie2 alignment length", default=50, type=int)
parser.add_argument("--min_accuracy", metavar="CONFIDENCE", help="Minimum threshold on Confidence score (percentage) to pass the reconstruction step", default=0.90, type=float)
parser.add_argument("--debug", help="Debug Mode", action='store_true', default=False) 
parser.add_argument("--presorted", help="The input BAM file is sorted and indexed with samtools. If set, MetaMLST skips this step", action='store_true') 
parser.add_argument("--quiet", help="Suppress text output", action='store_true') 
parser.add_argument("--legacy_samtools", help="Legacy mode (for samtools < 1.3)", action='store_true')
parser.add_argument("--version", help="Prints version informations", action='store_true')


parser.add_argument("--nloci", metavar="NLOCI", help="Do not discard samples where at least NLOCI (percent) are detected. This can lead to imperfect MLST typing", default=100, type=int)
parser.add_argument("--log", help="generate logfiles", action="store_true") 
parser.add_argument("-a", help="Write known sequences", action="store_true")

#parser.print_help()

args=parser.parse_args()


if args.version: print_version()

#PREPARE 
try:
	MetaMLSTDBconn = sqlite3.connect(args.d)
	MetaMLSTDBconn.row_factory = sqlite3.Row
	c = MetaMLSTDBconn.cursor()
except IOError: 
	metamlst_print("Failed to connect to the database: please check your database file!",'FAIL',bcolors.FAIL)
	sys.exit(1)




try:
	fil = open(args.BAMFILE,'r')
except IOError as e: 
	metamlst_print('Unable to open the BAM file for reading: please check the path is correct','FAIL',bcolors.FAIL)
	sys.exit(1)

cel = {}
sequenceBank = {}
ignoredReads = 0
totalReads = 0
glob_bamFile_sorted=False

fileName = (args.BAMFILE.split('/'))[-1].split('.')[0] if args.BAMFILE != '' else 'STDIN_'+str(int(time.time()))

if not os.path.isdir(args.o): os.mkdir(args.o)
workUnit = args.o	
 

try:
	child = subprocess.Popen(["samtools","view","-h","-"], stdout=subprocess.PIPE, stdin = fil)
except OSError as e:
	print "Error while executing samtools. Please check samtools is installed properly on your system!",e
	sys.exit(1)
 
for line in child.stdout:
	if(line[0] != '@'): 
		read = line.split('\t')
		readCode = read[0]

		species,gene,allele = read[2].split('_')
		
		score = (int)((read[11]).split(':')[2])
		xM = (int)((read[14]).split(':')[2])
		sequence = read[9]
		quality = read[10]
		
		if (args.filter and species in args.filter.split(',')) or not args.filter:
			if score >= args.minscore and len(sequence) >= args.min_read_len and xM <= args.max_xM:
				if species not in cel: cel[species] = {} #new species
			
				if gene not in cel[species]: #found a new gene
					cel[species][gene] = {} 
					sequenceBank[species+'_'+gene] = {}
					
				if allele not in cel[species][gene]: #found a new allele
					cel[species][gene][allele] = []
				
				cel[species][gene][allele].append(score)

				sequenceBank[species+'_'+gene][readCode]=len(sequence)

			else: ignoredReads = ignoredReads+1
			totalReads = totalReads +1

#COMPILES THE DATA STRUCTURE 			
for speciesKey,species in cel.items():

	for geneKey, geneInfo in species.items(): #geni 
	
		maxLen = max([len(x) for x in geneInfo.values()])

		for alleleCode,bowtieValues in geneInfo.items(): #alleli passata 2
			geneLen = len(bowtieValues)
			
			localScore = sum(item for item in bowtieValues)
			
			#print "GENELEN",geneKey,alleleCode,'a=',geneLen,'b=',maxLen, 'c=',localScore
			
			if geneLen != maxLen:
				localScore = localScore - (maxLen-geneLen)*args.penalty

			averageScore = float(localScore)/float(geneLen)

			cel[speciesKey][geneKey][alleleCode] = (localScore,geneLen,round(averageScore,1)) 
		

		

# cel ['haemophilus']['haemophilus_gene'][allele] = (a, b, c) --> a: summed score, b = geneLen, c = a/b

###################################### OUT FILE 1 (LOG)
if args.log:
	dfil = open(workUnit+'/'+fileName+'_'+str(int(time.time()))+'.out','w')
	dfil.write("SAMPLE:\t\t\t\t\t"+args.BAMFILE+'\r\n')
	dfil.write("PENALTY:\t\t\t\t"+repr(args.penalty)+'\r\n')
	dfil.write("MIN-THRESHOLD SCORE:\t\t\t\t"+repr(args.minscore)+'\r\n')
	
	dfil.write("TOTAL ALIGNED READS:\t\t\t\t"+repr(totalReads)+'\r\n')
	dfil.write(" - OF WHICH IGNORED:\t\t\t\t"+repr(ignoredReads)+' BAM READS\r\n\r\n------------------------------  RESULTS ------------------------------\r\n')

	for speciesKey,species in cel.items():
		for geneKey, geneInfo in species.items(): #geni
			for geneInfoKey,(score,geneLen,average) in sorted(geneInfo.items(), key= lambda x: x[1]): #alleli ordinati
				dfil.write('\t'.join(map(str,[speciesKey,geneKey,geneInfoKey,score,geneLen,average]))+'\r\n')
	dfil.close()	

 
if not args.quiet: print '\r\n'+bcolors.OKBLUE+('  '+fileName+'  ').center(80,'-')+bcolors.ENDC

for speciesKey,species in cel.items():
	
	#tVar is the dictionary of locus detected in the organism
	tVar = dict([(row['geneName'],0) for row in  c.execute("SELECT geneName FROM genes WHERE bacterium = ?",(speciesKey,))])
	
	#GENE PRESENCE 
	
	if len(tVar) < len(species.keys()):
		if not args.quiet: print 'Database is broken for' +speciesKey+bcolors.FAIL+'[ - EXITING - ]'.rjust(75,' ')+bcolors.ENDC
		sys.exit(0)
	
	for sk in species.keys():
		tVar[sk] = 1
	vals = sum([t for t in tVar.values()]) 
	
	if not args.quiet: 
		print ((bcolors.OKGREEN if (int((float(vals)/float(len(tVar)))*100) >= args.nloci) else bcolors.FAIL)+' '+speciesKey.ljust(18,' ')+bcolors.ENDC)+' Detected Loci: '+', '.join([bcolors.OKGREEN + sk + bcolors.ENDC for (sk,v) in sorted(tVar.items(), key = lambda x:x[0]) if v == 1])
		if len([ta for ta in tVar.items() if v == 0]) > 0: print (' '*20)+'Missing Loci : '+', '.join([bcolors.FAIL + sk + bcolors.ENDC for (sk,v) in sorted(tVar.items(), key = lambda x:x[0]) if v == 0])
		print ""
	
	if int((float(vals)/float(len(tVar)))*100) >= args.nloci:
	
		if not args.quiet: metamlst_print("Closest allele identification",'...',bcolors.HEADER)

		if not args.quiet: print "\r\n  "+"Locus".ljust(7)+"Avg. Coverage".rjust(15)+"Score".rjust(7)+"Hits".rjust(6)+" Reference Allele(s)".ljust(36)
		sys.stdout.flush()
		
		for geneKey, geneInfo in sorted(species.items(),key=lambda x:x[0]): #loci
		
			minValue = max([avg for (val,leng,avg) in geneInfo.values()])
			
			aElements = {}
			for k,(val,leng,avg) in geneInfo.items():
				if avg == minValue:
					aElements[k]=(val,leng,avg)

			closeAllelesList = ','.join([str(a) for a in sorted(aElements.keys(), key=lambda x: int(x))][:5])+('... ('+str(len(aElements))+' more)'  if len(aElements) > 5 else '')
		
			sequenceKey = speciesKey+'_'+geneKey
			c.execute("SELECT LENGTH(sequence) as L FROM alleles WHERE bacterium = ? AND gene = ? ORDER BY L DESC LIMIT 1", (speciesKey,geneKey))
			genL = c.fetchone()['L']

			coverage = sum([x for x in sequenceBank[sequenceKey].values()])

			if not args.quiet: print "  "+bcolors.WARNING+geneKey.ljust(7)+bcolors.ENDC+str(round(float(coverage)/float(genL),2)).rjust(15)+bcolors.ENDC+bcolors.HEADER+str(minValue).rjust(7)+str(aElements.itervalues().next()[1]).rjust(6)+bcolors.ENDC+bcolors.OKBLUE,closeAllelesList.ljust(36)+bcolors.ENDC
		if not args.quiet: print ""
		sys.stdout.flush()		
		
		sampleSequence = '' #for organism
		
		if not args.quiet: metamlst_print("Building Consensus Sequences",'...',bcolors.HEADER)
		

		if not args.presorted and not glob_bamFile_sorted:
			
			sort_index(args.BAMFILE,legacy=args.legacy_samtools)
			glob_bamFile_sorted=True

		l = [sorted([(speciesKey+'_'+g1+'_'+k,db_getUnalSequence(MetaMLSTDBconn,speciesKey,g1,k)) for k,(val,leng,avg) in g2.items() if avg == max([avg1 for (val1,leng1,avg1) in g2.values()])],key=lambda x: int(x[0].split('_')[2]))[0] for g1,g2 in species.items()] 
		consenSeq = buildConsensus(args.BAMFILE, dict(l),args.minscore,args.max_xM,args.debug,legacy=args.legacy_samtools)
		
		
		
		newProfile = 0
		finWrite = 1
		if not args.quiet: print "\r\n  "+"Locus".ljust(7)+"Ref.".ljust(7)+"Length".rjust(7)+"Ns".rjust(7)+"SNPs".rjust(7)+"Confidence".rjust(15)+"Notes".rjust(10)
		for l in sorted(consenSeq, key= lambda x: x.id):
			holes = str(l.description.split('_')[0].split('::')[1])
			snps = int(l.description.split('_')[1].split('::')[1])
			leng = str(len(l.seq))
			leng_ns = str(round(1-float(holes)/float(leng),4)*100)+' %'
			l.seqLen = len(l.seq)
			
			# If there's a gene with low accuracy, the whole organism is discarded for this sample
			if (1-float(holes)/float(leng)) <= args.min_accuracy: finWrite = 0
			
			if snps > 0:
				seqFind = sequenceFind(MetaMLSTDBconn,speciesKey,l.seq)
				if seqFind: newAllele = seqFind
				else:
					newAllele = 'NEW'
					newProfile = 1
					
			else:
				newAllele = '--'
				if not args.a: l.seq = ''
			
			if not args.quiet: print "  "+bcolors.WARNING+(l.id.split('_')[1]).ljust(7)+bcolors.ENDC+(l.id.split('_')[2]).ljust(7)+leng.rjust(7)+holes.rjust(7)+str(snps).rjust(7)+leng_ns.rjust(15)+newAllele.rjust(10)
			
		
		if not args.quiet: print ''
		 
		
		if finWrite:
			
			if not args.quiet: metamlst_print("Reconstruction Successful",'WRITE',bcolors.OKGREEN)
			profil = open(workUnit+'/'+fileName+'.nfo','a')	 
			profil.write(speciesKey+'\t'+fileName+'\t'+"\t".join([recd.id+"::"+str(recd.seq)+'::'+str(round(1-float(recd.description.split('_')[0].split('::')[1])/float(recd.seqLen),4)*100)+'::'+str(round(float(recd.description.split('_')[1].split('::')[1])/float(recd.seqLen),4)*100) for recd in consenSeq])+'\r\n')
			
			profil.close()
		else:
			if not args.quiet: metamlst_print("Accuracy lower than "+str(round(args.min_accuracy*100,2))+'%','SKIP',bcolors.FAIL) 
		
	del cel[speciesKey]
	gc.collect()
			
 
MetaMLSTDBconn.close() 

if len(cel) and not args.quiet: print '\033[92m'+'[ - Completed - ]'.rjust(80,' ')+'\033[0m'

