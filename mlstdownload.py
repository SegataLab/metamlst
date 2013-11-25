#python

from Bio import SeqIO
import sqlite3
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import argparse
import urllib, urllib2, cookielib, re, os
from HTMLParser import HTMLParser
# create a subclass and override the handler methods
class MyHTMLParser(HTMLParser):
    def handle_data(self, data):
        if data.strip() != '': loci.append(data.strip())
htmlparser = MyHTMLParser()


parser = argparse.ArgumentParser()
parser.add_argument("bacteria", help="Bacteria subdomain in mlst website, COMMA SEPARATED")
parser.add_argument("-d","--download", help="Download Loci", action="store_true") 
args=parser.parse_args()

os.remove('bsb.db')
conn = sqlite3.connect('bsb.db')
c = conn.cursor()
c.execute("CREATE TABLE IF NOT EXISTS organisms (bacteriumName varchar(255) PRIMARY KEY)")
c.execute("CREATE TABLE IF NOT EXISTS genes (geneName varchar(255), bacterium VARCHAR(255), PRIMARY KEY(geneName,bacterium))")
c.execute("CREATE TABLE IF NOT EXISTS alleles (recID INTEGER PRIMARY KEY AUTOINCREMENT,alleleName varchar(255), gene VARCHAR(255), sequence TEXT, alleleVariant INT)")
c.execute("CREATE TABLE IF NOT EXISTS profiles (redID INTEGER PRIMARY KEY AUTOINCREMENT, profileCode INTEGER, species VARCHAR(255))")

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

		#FASTA FILE (bowtie index)
		out_file = open(organism+".faa","w")
		for locus in loci:
			url_2 = 'http://'+organism+'.mlst.net/sql/fasta.asp?allele='+locus
			req = urllib2.Request(url_2)
			rsp = urllib2.urlopen(req)
			content = rsp.read()
			cs = re.split('<textarea[^>]*>',content)
			cs = re.split('</textarea>',cs[1])
			out_file.write(cs[0].replace('&gt;','>'+organism+'_'))

		out_file.close()
		seqList = []
		#DATABASE FILE
		
		for seq_record in SeqIO.parse(organism+".faa", "fasta"):
			seqList.append(SeqRecord(Seq(str(seq_record.seq), IUPAC.unambiguous_dna), id = seq_record.id, description=''))  
			print str(seq_record.id)
			c.execute("INSERT INTO alleles (alleleName, gene,sequence, alleleVariant) VALUES (?,?,?,?)", (str(seq_record.id), re.findall('_[a-zA-Z_]*',str(seq_record.id))[0].replace('_',''), str(seq_record.seq),re.findall('[0-9]*$',str(seq_record.id))[0]))
		
		
		#PRIFILE FILE 
		url_2 = 'http://'+organism+'.mlst.net/sql/st_comma.asp'
		req = urllib2.Request(url_2)
		rsp = urllib2.urlopen(req)
		content = rsp.read()
		content.split('\r\n')
		print content
		
		#c.executemany("INSERT INTO alleles (alleleName, gene,sequence, alleleVariant) VALUES (?,?,?,?)", [(str(seq_record.id), re.findall('_[a-zA-Z_]*',str(seq_record.id))[0].replace('_',''), str(seq_record.seq),re.findall('[0-9]*$',str(seq_record.id))[0]) for seq_record in seqList]) 
 


conn.commit()
conn.close() 