#!/usr/bin/env python

import json, math,sqlite3,itertools,subprocess,time,gc,argparse, re, os,sys,pkgutil,importlib

from StringIO import StringIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


__author__ = 'Moreno Zolfo (moreno.zolfo@unitn.it)'
__version__ = '1.1'
__date__ = '26 April 2017'


def print_version():
	print "Version:\t"+__version__
	print "Author:\t\t"+__author__
	print "Reference:\t"+'MetaMLST: multi-locus strain-level bacterial typing from metagenomic samples\n\t\tNucleic Acids Research, 2016\n\t\tDOI: 10.1093/nar/gkw837'
	sys.exit(0)

def metamlst_print(mesg,label,type,reline=False,newLine=False):
	opening = "\r" if reline else ''
	ending = "\r\n" if not reline or newLine else ''

	if len(mesg) < 65:
	
		sys.stdout.write(opening+mesg.ljust(66)+(type+'[ - '+label.center(5)+' - ]'+bcolors.ENDC).ljust(14)+ending)
	else: 
		c=0
		wds = []
		lines=[]
		for word in mesg.split(' '):

				if c + len(word)+2 > 65:
					print ' '.join(wds)
					c=0
					wds=[word]
					continue
				c = c+len(word)+2
				wds.append(word)
		sys.stdout.write(opening+(' '.join(wds)).ljust(66)+(type+'[ - '+label.center(5)+' - ]'+bcolors.ENDC).ljust(14)+ending)

	sys.stdout.flush()

def metamlst_newline():
	sys.stdout.write('\r\n')

def dump_db_to_fasta(conn,path):
	conn.row_factory = sqlite3.Row
	cursor = conn.cursor() 
	
	#print '  COLLECTING DATA'.ljust(60),
	metamlst_print("COLLECTING DATA...",'...',bcolors.ENDC)
	sys.stdout.flush()
	
	strr=[SeqRecord(Seq(row['sequence'],IUPAC.unambiguous_dna),id=row['bacterium']+'_'+(row['gene']+'_'+str(row['alleleVariant'])),description='') for row in cursor.execute("SELECT bacterium,gene,alleleVariant,sequence FROM alleles WHERE sequence <> ''")]	
	SeqIO.write(strr,path,'fasta')
	return len(strr)

def alleleInProfile(conn,bacterium,gene,allele,profile):
	e = conn.cursor()
	e.execute("SELECT 1 FROM profiles,alleles WHERE alleleCode = alleles.recID AND profileCode = ? AND alleleVariant = ? AND gene = ? AND profiles.bacterium = ?",(profile,allele,gene,bacterium))
	return (len(e.fetchall()) > 0)
	
def sequenceExists(conn,bacterium,sequence):
	
	e = conn.cursor()
	e.execute("SELECT 1 FROM alleles WHERE sequence = ? AND bacterium = ?",(str(sequence),bacterium))
	return (len(e.fetchall()) > 0)

def alleleInManyProfile(conn,bacterium,gene,allele,profileList):
	
	e = conn.cursor()
	profileString = ','.join([str(x) for x in profileList])
	e.execute("SELECT 1 FROM profiles,alleles WHERE alleleCode = alleles.recID AND profileCode IN ("+profileString+")  AND alleleVariant = ? AND gene = ? AND profiles.bacterium = ?",(allele,gene,bacterium))
	return (len(e.fetchall()) > 0)
	
def db_getSequence(conn,bacterium,gene,allele):
	e = conn.cursor()
	e.execute("SELECT alignedSequence FROM alleles WHERE bacterium = ? AND gene = ? AND alleleVariant = ?",(bacterium,gene,allele))
	return (e.fetchone()['alignedSequence'])

def db_getUnalSequence(conn,bacterium,gene,allele):
	e = conn.cursor()
	unal = e.execute("SELECT sequence FROM alleles WHERE bacterium = ? AND gene = ? AND alleleVariant = ?",(bacterium,gene,allele))
	unalobj = unal.fetchone()
	if unalobj is not None:
		return unalobj['sequence']
	else:
		metamlst_print(' > '+bacterium+'_'+gene+'_'+allele+" was not found in the database!",'!!!',bcolors.WARNING)
		return None

def sequenceFind(conn,bacterium,sequence):
	
	e = conn.cursor()
	e.execute("SELECT gene,alleleVariant FROM alleles WHERE sequence = ? AND bacterium = ?",(str(sequence),bacterium))
	res = e.fetchone()
	
	if res: return res['gene']
	else: return 0

def defineProfile(conn,geneList):
	
	recs=[]
	for allele in geneList:
		e = conn.cursor()
		e.execute("SELECT recID FROM alleles WHERE bacterium||'_'||gene||'_'||alleleVariant = ?",(allele,))
		result = e.fetchone()
		
		if result:
			recs.append(str(result['recID']))

	return [(row['profileCode'],int((float(row['T'])/float(len(recs)))*100)) for row in e.execute("SELECT profileCode, COUNT(*) as T FROM profiles WHERE alleleCode IN ("+','.join(recs)+") GROUP BY profileCode HAVING T = (SELECT COUNT(*)  FROM profiles WHERE alleleCode IN ("+','.join(recs)+") GROUP BY profileCode ORDER BY COUNT(*) DESC LIMIT 1) ORDER BY T DESC")] if result else [(0,0)]

def sequenceLocate(conn,bacterium,sequence):
	
	e = conn.cursor()
	e.execute("SELECT alleleVariant FROM alleles WHERE sequence = ? AND bacterium = ?",(str(sequence),bacterium))
	return str(e.fetchone()['alleleVariant'])

def sequencesGetAll(conn,bacterium,gene):
	
	e = conn.cursor()
	e.execute("SELECT sequence,alleleVariant FROM alleles WHERE gene = ? AND bacterium = ?",(gene,bacterium))
	return dict((x['alleleVariant'],x['sequence']) for x in e.fetchall())

def stringDiff(s1,s2):
	c =0 
	for a,b in zip(s1,s2):
		if a!=b: c+=1
	return c 


def sort_index(bamFile,legacy=False):
	if legacy:

		subprocess.call(['samtools','sort',bamFile,bamFile+'.s'])
		os.rename(bamFile+'.s.bam',bamFile)

	else:
		subprocess.call(['samtools','sort',bamFile,'-o',bamFile+'.sorted'])
		os.rename(bamFile+'.sorted',bamFile)
	
	subprocess.call(['samtools','index',bamFile])

def buildConsensus(bamFile,chromosomeList,filterScore,max_xM,debugMode,legacy=False):

	if legacy: metamlst_print('Samtools < 1.x legacy mode enabled','...',bcolors.HEADER)

 	if pkgutil.find_loader('pysam') is not None:
 		import pysam
	else:
		return buildConsensus_legacy(bamFile,chromosomeList,filterScore,max_xM,debugMode,legacy_samtools=legacy)

	if not os.path.isfile(bamFile+'.bai'): pysam.index(bamFile)
	seqRec = []
 	chromosomes = {}

	samfile = pysam.AlignmentFile(bamFile, "rb")

	chromosomesLen = dict((r,l) for r,l in zip(samfile.references,samfile.lengths))

	for chre in chromosomeList.keys():
		
		chromosomes[chre] = {}

		for pileupcolumn in samfile.pileup(chre,stepper='nofilter'):
			ldict = {'A':0,'T':0,'C':0,'G':0}
			#otherBases = {}			

			for pileupread in pileupcolumn.pileups:
				if not pileupread.is_del and not pileupread.is_refskip:
					base = pileupread.alignment.query_sequence[pileupread.query_position]
					
					if int(pileupread.alignment.get_tag('AS')) >= filterScore and int(pileupread.alignment.get_tag('XM')) <= max_xM: 
						
						if base in ldict.keys(): ldict[base]+=1
					#else otherBases.append(pileupread.alignment.query_sequence[pileupread.query_position])

			#print chre,pileupcolumn.pos+1,ldict
			chromosomes[chre][pileupcolumn.pos+1] = max(ldict, key=ldict.get)

	
	for chromo,nucleots in chromosomes.items(): 
		sys.stdout.flush()

		sequen = ""
		lastSet = 0

		for key,nucleotide in sorted(nucleots.items(), key=lambda x : x[0]):
			if lastSet+1 != key:
				sequen = sequen + "N" * (key-(lastSet+1))
			lastSet = key
			sequen = sequen + nucleotide
			
		sequen+=((chromosomesLen[chromo]-len(sequen))*"N")
		
		rSequen = list(sequen)
		dbSequen = chromosomeList[chromo] 
		
		i=0
		cIndex=0
		SNPs=0 
		
		
		for chr in rSequen:
			if chr == 'N':
				rSequen[i] = dbSequen[i].lower()
				cIndex+=1
				
			elif rSequen[i] != dbSequen[i]: 
				rSequen[i] = rSequen[i]
				SNPs+=1
			i+=1
		sequen = ''.join(rSequen)
		
			
		seqRec.append(SeqRecord(Seq(sequen,IUPAC.unambiguous_dna),id=chromo, description = 'CI::'+str(cIndex)+'_SP::'+str(SNPs)))
	print '\r',


	samfile.close()
	#seqRec.append(SeqRecord(Seq(sequen,IUPAC.unambiguous_dna),id=chromo, description = 'CI::'+str(cIndex)+'_SP::'+str(SNPs)))
	
	return seqRec
def buildConsensus_legacy(bamFile,chromosomeList,filterScore,max_xM,debugMode,legacy_samtools=False):

	chromosomes = {}
	chromosomesLen = {}
	
	metamlst_print('Pysam not detected: using legacy mode','...',bcolors.HEADER)
	strr=''
	filterSpec = chromosomeList.keys()[0].split('_')[0].strip()
	
	fili = open(bamFile,'r') if os.path.isfile(bamFile) else sys.stdin
	devnull = open('/dev/null', 'w')
	child = subprocess.Popen("samtools view -h - | grep "+filterSpec,shell=True, stdout=subprocess.PIPE, stdin = fili)
	
	for line in child.stdout:
	
		pox = 0 
		line = line.strip()
		
		if line=='': continue
		
		for chromoKey in chromosomeList.keys():
			if  (line.startswith('@SQ') and line.split('\t')[1].split(':')[1].strip() == chromoKey) or (not line.startswith('@') and line.split('\t')[2] == chromoKey):
					pox = 1 
					break 
				
		if not pox: 
			continue
		
		if line.startswith("@SQ"): 						#chromosomes data
			chromosomesLen[line.split('\t')[1].replace('SN:','')] = int(line.split('\t')[2].replace('LN:',''))
			
			strr+=line+'\r\n'
		elif line.startswith("@"): strr+=line+'\r\n'	#other data
		else: 											#align data
			tline = line.split('\t')
			scoreOfRead = (int)((tline[11]).split(':')[2])
			xM = (int)((tline[14]).split(':')[2])
			
			if scoreOfRead >= filterScore and xM <= max_xM:
				tline[1] = str(int(tline[1]) & 0b1111011111111)
				strr+='\t'.join(tline)+'\r\n' 
	
	
	devnull = open('/dev/null', 'w')
	
	if legacy_samtools: child = subprocess.Popen("samtools view -bS - | samtools sort -o - - | samtools mpileup - ",shell=True, stdout=subprocess.PIPE, stdin = subprocess.PIPE,stderr = devnull)
	else: child = subprocess.Popen("samtools view -bS - | samtools sort -o - - | samtools mpileup - ",shell=True, stdout=subprocess.PIPE, stdin = subprocess.PIPE,stderr = devnull)

	out = StringIO(child.communicate(strr)[0]) 
	
	if debugMode: print out.getvalue() 
	
	for line in out:  
		chromosome = line.split('\t')[0]
		nucleotide = int(line.split('\t')[1])  
		
		if line.split('\t')[4] == '': continue
		
		escape = 0		
		escapeLen = 0	
		
		if chromosome not in chromosomes: chromosomes[chromosome] = {}
		ldict = {'A':0, 'T':0, 'C':0,'G':0}
		for chr in line.split('\t')[4]:
			
			
			if chr == '+' or chr == '-': 
				escape = 1
				continue
				
			if escape: 
				if chr.isdigit(): 
					escapeLen = int(chr)+1
				escape = 0 
				
			if escapeLen > 0: 
				escapeLen -=1
				continue
			
			elif chr.upper() in ldict: 
				ldict[chr.upper()] += 1
			
		chromosomes[chromosome][nucleotide] = max(ldict, key=ldict.get)
	
	del out
	seqRec=[]

	for chromo,nucleots in chromosomes.items(): 
		sys.stdout.flush()
		
		sequen = ""
		lastSet = 0

		for key,nucleotide in sorted(nucleots.items(), key=lambda x : x[0]):
			if lastSet+1 != key:
				sequen = sequen + "N" * (key-(lastSet+1))
			lastSet = key
			sequen = sequen + nucleotide
			
		sequen+=((chromosomesLen[chromo]-len(sequen))*"N")
		
		rSequen = list(sequen)
		dbSequen = chromosomeList[chromo] 
		
		i=0
		cIndex=0
		SNPs=0 
		
		
		for chr in rSequen:
			if chr == 'N':
				rSequen[i] = dbSequen[i].lower()
				cIndex+=1
				
			elif rSequen[i] != dbSequen[i]: 
				rSequen[i] = rSequen[i]
				SNPs+=1
			i+=1
		sequen = ''.join(rSequen)
		
			
		seqRec.append(SeqRecord(Seq(sequen,IUPAC.unambiguous_dna),id=chromo, description = 'CI::'+str(cIndex)+'_SP::'+str(SNPs)))
	print '\r',
	return seqRec

MLST_KEYWORDS = ['clonal_complex','species','mlst_clade']
	
class bcolors:
	HEADER = '\033[95m'
	OKBLUE = '\033[94m'
	OKGREEN = '\033[92m'
	WARNING = '\033[93m'
	FAIL = '\033[91m'
	ENDC = '\033[0m'
	OKGREEN2 = '\033[42m\033[30m'
	RED = '\033[1;91m'
	CYAN = '\033[0;37m'

def db_getOrganisms(conn,bacterium=None):
	e = conn.cursor()
	t= dict((elem['organismkey'],(elem['label']) if elem['label'] is not None else '('+elem['organismkey']+')') for elem in e.execute("SELECT * FROM organisms"))
	if bacterium: return t[bacterium] 
	else: return  t
	
class metaMLST_db:
	
	cursor = None
	conn = None

	def __init__(self,dbPath):
		
		try:
			self.conn = sqlite3.connect(dbPath)
			self.conn.row_factory = sqlite3.Row
			self.cursor = self.conn.cursor()
		except IOError: 
			print "IOError: unable to access "+args.database+"!"
	
	def getOrganisms(self,bacterium=None):
		listAlleles = []

		t= dict((elem['organismkey'],(elem['label']) if elem['label'] is not None else '('+elem['organismkey']+')') for elem in self.cursor.execute("SELECT * FROM organisms"))
		if bacterium: return t[bacterium] 
		else: return  t
	

	def getAlleles(self,profile):
		listAlleles = []
		for row in self.cursor.execute("SELECT bacterium,gene,alleleVariant,sequence FROM alleles WHERE sequence <> '' AND bacterium = ?",(profile,)):
			listAlleles.append(SeqRecord(Seq(row['sequence'],IUPAC.unambiguous_dna),id=row['bacterium']+'_'+(row['gene']+'_'+str(row['alleleVariant'])),description=''))
		return listAlleles
		
	def getGene(self,profile,geneName):
		listAlleles = []
		for row in self.cursor.execute("SELECT bacterium,gene,alleleVariant,sequence FROM alleles WHERE sequence <> '' AND bacterium = ? AND gene = ?",(profile,geneName)):
			listAlleles.append(SeqRecord(Seq(row['sequence'],IUPAC.unambiguous_dna),id=row['bacterium']+'_'+(row['gene']+'_'+str(row['alleleVariant'])),description=''))
		return listAlleles
		
	def getGeneNames(self,profile):
		lister = []
		for row in self.cursor.execute("SELECT geneName FROM genes WHERE bacterium = ?",(profile,)):
			lister.append(row['geneName'])
		return lister
		

	def defineProfile(self,geneList):
	
		recs=[]
		result = None
		for allele in geneList:
			self.cursor.execute("SELECT recID FROM alleles WHERE bacterium||'_'||gene||'_'||alleleVariant = ?",(allele,))
			result = self.cursor.fetchone()
			if result: recs.append(str(result['recID']))
		
		return [(row['profileCode'],int((float(row['T'])/float(len(recs)))*100)) for row in self.cursor.execute("SELECT profileCode, COUNT(*) as T FROM profiles WHERE alleleCode IN ("+','.join(recs)+") GROUP BY profileCode HAVING T = (SELECT COUNT(*)  FROM profiles WHERE alleleCode IN ("+','.join(recs)+") GROUP BY profileCode ORDER BY COUNT(*) DESC LIMIT 1) ORDER BY T DESC")] if result else [(0,0)]