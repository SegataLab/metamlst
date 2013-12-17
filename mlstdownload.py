#python

from Bio import SeqIO
import sqlite3
from StringIO import StringIO
from Bio.Seq import Seq
from Bio.Align.Applications import MuscleCommandline
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import argparse
import urllib, urllib2, cookielib, re, os, subprocess, sys
from HTMLParser import HTMLParser
# create a subclass and override the handler methods
class MyHTMLParser(HTMLParser):
    def handle_data(self, data):
        if data.strip() != '': loci.append(data.strip())
htmlparser = MyHTMLParser()

class bcolors:
	HEADER = '\033[95m'
	OKBLUE = '\033[94m'
	OKGREEN = '\033[92m'
	WARNING = '\033[93m'
	FAIL = '\033[91m'
	ENDC = '\033[0m'
	OKGREEN2 = '\033[42m\033[30m'

parser = argparse.ArgumentParser()
parser.add_argument("bacteria", help="Bacteria subdomain in mlst website, COMMA SEPARATED")
parser.add_argument("-d","--download", help="Download Loci", action="store_true") 
parser.add_argument("-stf","--stfile", help="GET STs from COMMA SEPARATED FILE") 
args=parser.parse_args()

#os.remove('bsb.db')
conn = sqlite3.connect('bsb.db')
conn.row_factory = sqlite3.Row
c = conn.cursor()
c.execute("CREATE TABLE IF NOT EXISTS organisms (bacteriumName varchar(255) PRIMARY KEY)")
c.execute("CREATE TABLE IF NOT EXISTS genes (geneName varchar(255), bacterium VARCHAR(255), PRIMARY KEY(geneName,bacterium))")
c.execute("CREATE TABLE IF NOT EXISTS alleles (recID INTEGER PRIMARY KEY AUTOINCREMENT,bacterium varchar(255), gene VARCHAR(255), sequence TEXT, alignedSequence TEXT, alleleVariant INT)")
c.execute("CREATE TABLE IF NOT EXISTS profiles (recID INTEGER PRIMARY KEY AUTOINCREMENT, profileCode INTEGER, bacterium VARCHAR(255), alleleCode INTEGER)")

#CONN TO MLST
print args.bacteria.split(',')
for organism in args.bacteria.split(','):
	url_2 = 'http://'+organism+'.mlst.net/sql/download_alleles.asp'
	req = urllib2.Request(url_2)
	rsp = urllib2.urlopen(req)
	content = rsp.read()
	cs = re.split('<select name="allele"[^>]*>',content)
	cs = re.split('</select>',cs[1])
	loci = []
	htmlparser.feed(cs[0]);

	print "Loci for "+organism+" are:\t",loci 

	c.execute("INSERT INTO organisms (bacteriumName) VALUES (?)",(organism,))
	c.executemany("INSERT INTO genes (geneNAme, bacterium) VALUES (?,?)", [(locus,organism) for locus in loci])
 
	if args.download:
		seqToAlign = {}
	
		#FASTA FILE (bowtie index)
		out_file = open(organism+".faa","w")
		for locus in loci:
			print bcolors.OKBLUE+locus+bcolors.ENDC
			print(('\tDownloading genes/alleles:').ljust(60)),
			url_2 = 'http://'+organism+'.mlst.net/sql/fasta.asp?allele='+locus
			req = urllib2.Request(url_2)
			rsp = urllib2.urlopen(req)
			content = rsp.read()
			cs = re.split('<textarea[^>]*>',content)
			cs = re.split('</textarea>',cs[1])
			out_file.write(cs[0].replace('&gt;','>'+organism+'_')) #write for bowtie index .fasta file
			print(bcolors.OKGREEN+'[ - Done - ]'+bcolors.ENDC)
			#from StringIO import StringIO
			cString = StringIO(cs[0].replace('&gt;','>'+organism+'_'))
			
			recList = SeqIO.parse(cString, "fasta")
			
			#ALIGNMENT
			print(('\tAligning Sequences:').ljust(60)),
			sys.stdout.flush()
			muscle_cline = MuscleCommandline("tools/muscle")
			alignTable={}
			stdout, stderr = muscle_cline(stdin=cString.getvalue())
			for sequence in SeqIO.parse(StringIO(stdout), "fasta"):
				alignTable[sequence.id] = str(sequence.seq)
			print(bcolors.OKGREEN+'[ - Done - ]'+bcolors.ENDC)
			
			#INSERTION
			insertionList = []
			for seq_record in recList:
				#seqList.append(SeqRecord(Seq(str(seq_record.seq), IUPAC.unambiguous_dna), id = seq_record.id, description=''))  
				alleleName = str(seq_record.id)
				gene = re.sub('^_','',re.findall('_[a-zA-Z_]*',str(seq_record.id))[0])
				sequence =  str(seq_record.seq)
				alleleVariant = re.findall('[0-9]*$',str(seq_record.id))[0]
				insertionList.append((gene,organism,sequence,alignTable[alleleName],alleleVariant))
				
			print(('\tWriting to DB:').ljust(60)),
			c.executemany("INSERT INTO alleles (gene, bacterium,sequence,alignedSequence, alleleVariant) VALUES (?,?,?,?,?)", insertionList)
			print(bcolors.OKGREEN+'[ - Done - ]'+bcolors.ENDC)
			

		out_file.close()
		seqList = []
		
		#PROFILE FILE 
		print(('Parsing MLST profiles:').ljust(67)),
		if args.stfile: 
			fil = open(args.stfile,'r')
			content = fil.read()
			content = content.split('\n')
		else:
			url_2 = 'http://'+organism+'.mlst.net/sql/st_comma.asp'
			req = urllib2.Request(url_2)
			rsp = urllib2.urlopen(req)
			content = rsp.read()
			content = content.split('<br>')
		
		for profile in content:
			profile = profile.split(',')
			
			if len(profile) > 1:
				profileCode = profile[0]
				cI = 0
				insertion = []
				for row in c.execute("SELECT geneName FROM genes WHERE bacterium = ?",(organism,)):
					cI = cI+1
					insertion.append((profileCode,organism,organism,str(row['geneName']),profile[cI]))
					
				c.executemany("INSERT INTO profiles (profileCode, bacterium, alleleCode) VALUES (?,?,(SELECT recID FROM alleles WHERE bacterium = ? AND gene = ? AND alleleVariant = ?))", insertion)
						
		print(bcolors.OKGREEN+'[ - Done - ]'+bcolors.ENDC)	

conn.commit()
conn.close() 