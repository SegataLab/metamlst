#!/usr/bin/env python


from Bio.Align.Applications import MuscleCommandline
from Bio import SeqIO
import sqlite3
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import argparse
import metaMLST_functions 

from metaMLST_functions import * 




parser = argparse.ArgumentParser()

parser.add_argument("database", help="database file")
parser.add_argument("--cli_correct_brutal", help="corrects (and writes!) on the db (remove aberrant organisms if --cli_correct can't fix them)", action="store_true") 
parser.add_argument("--cli_correct", help="corrects (and writes!) on the db (remove aberrant alleles with .9 threshold)", action="store_true") 
parser.add_argument("--cli_correct_except", help="with cli_correct, extends the correction to a specific organism regardless of .9 threshold") 
parser.add_argument("--cli", help="checks sequences in database", action="store_true") 
parser.add_argument("--probe_gene", help="Visualize all about a gene Must be in the format: species_gene_") 
parser.add_argument("--remove_allele", help="removes an allele from the database. Must be in the format: species_gene_allelenumber") 
parser.add_argument("--remove_gene", help="removes a gene and all its alleles from the database. Must be in the format: species_gene") 
parser.add_argument("--log", help="generates a logfile", action="store_true") 

args = parser.parse_args()

try:
  conn = sqlite3.connect(args.database)
  conn.row_factory = sqlite3.Row
  cursor = conn.cursor() 
except IOError:
  print "IOError: unable to access "+args.database+"!"

conn.row_factory = sqlite3.Row
cursor = conn.cursor() 

if args.log: logFile=open('log.log','w')

if args.probe_gene: 
  (organism,gene) = args.probe_gene.split('_')
  print "ID\tGENE\tALLELE\tSEQ"
  for elem in cursor.execute("SELECT * FROM alleles WHERE bacterium = ? AND gene = ?",(organism,gene)):
    print elem['recID'],'\t',elem['gene'],'\t',elem['alleleVariant'],elem['sequence']
  

if args.remove_allele:
  (organism,gene,allele) = args.remove_allele.split('_')
  if len([row for row in cursor.execute("SELECT * FROM alleles WHERE bacterium = ? AND gene = ? AND alleleVariant = ?",(organism,gene,allele))]) > 0:
    cursor.execute("DELETE FROM alleles WHERE bacterium = ? AND gene = ? AND alleleVariant = ?",(organism,gene,allele))

if args.remove_gene:
  (organism,gene) = args.remove_gene.split('_')
  if len([row for row in cursor.execute("SELECT 1 FROM alleles WHERE bacterium = ? AND gene = ?",(organism,gene))]) > 0: cursor.execute("DELETE FROM alleles WHERE bacterium = ? AND gene = ?",(organism,gene))
  if len([row for row in cursor.execute("SELECT 1 FROM genes WHERE bacterium = ? AND geneName = ?",(organism,gene))]) > 0: cursor.execute("DELETE FROM genes WHERE bacterium = ? AND geneName = ?",(organism,gene))
  
if args.cli:
  base = {}
  for x in cursor.execute("SELECT * FROM alleles"):
      organism = x['bacterium']
      gene = x['gene'] 
      if organism not in base: base[organism] = {}
      if gene not in base[organism]: base[organism][gene] = []
      base[organism][gene].append((x['alleleVariant'],x['sequence']))
  
  for bacteriumName,bacteriumRec in base.items():
      print bcolors.OKBLUE+bacteriumName+bcolors.ENDC
      for geneName,geneList in bacteriumRec.items():
	 if len(set([len(x) for (a,x) in geneList])) > 1:
	  print '\t',bcolors.WARNING+geneName+bcolors.ENDC
	  cntl = {}
	  for (variant,seq) in geneList:
	      if len(seq) not in cntl: cntl[len(seq)] = 0
	      cntl[len(seq)] +=1
	  
	  total_alleles=sum(cntl.values())
	  len_of_max = max(cntl,key=cntl.get)
	  per_of_max= float(cntl[len_of_max]) / float(total_alleles)
	  
	  # print '\t\t',(str(len_of_max)+' bps').ljust(10)+(str(round(float(cntl[len_of_max]) / float(total_alleles),4)*100).rjust(5)+'% *').ljust(12)+str(cntl[len_of_max]).rjust(4)+' / '+str(total_alleles).rjust(4)
	  for removableLen, removableCount in sorted(cntl.items(),key=lambda x:x[1],reverse=True):
	    if removableLen == len_of_max:
	      bcl= bcolors.OKGREEN 
	      remFlag = ''
	    else:
	      bcl = ''
	      remFlag = bcolors.FAIL+'REMOVE'+bcolors.ENDC if per_of_max >= 0.9 else bcolors.WARNING+'CHECK'+bcolors.ENDC
		 
	    print '\t\t'+bcl+(str(removableLen)+' bps').ljust(10)+(str(round(float(removableCount) / float(total_alleles),4)*100).rjust(5)+'%').ljust(8)+str(removableCount).rjust(4)+' / '+str(total_alleles).rjust(4)+' '+remFlag.rjust(8)+bcolors.ENDC
	  
	  #end of gene
	  if args.cli_correct and (per_of_max >= 0.9 or bacteriumName == args.cli_correct_except): 
	      #print [row for row in cursor.execute("SELECT gene,alleleVariant,LENGTH(sequence) FROM alleles WHERE LENGTH(sequence) <> ? AND bacterium = ? and gene = ?", (len_of_max,bacteriumName,geneName))]
	      cursor.execute("DELETE FROM alleles WHERE LENGTH(sequence) <> ? AND bacterium = ? and gene = ?", (len_of_max,bacteriumName,geneName))
	      print ''.ljust(60),bcolors.OKGREEN+'[ - CORRECTED - ]'+bcolors.ENDC
	      if args.log: logFile.write(bacteriumName+'_'+geneName+':\t fixed\n')
	    #print str(cntl)+'\t',per_of_max, 
	  elif args.cli_correct_brutal and per_of_max < 0.9:
	      cursor.execute("DELETE FROM organisms WHERE organismkey = ?", (bacteriumName,))
	      cursor.execute("DELETE FROM genes WHERE bacterium = ?", (bacteriumName,))
	      cursor.execute("DELETE FROM alleles WHERE bacterium = ?", (bacteriumName,))
	      cursor.execute("DELETE FROM profiles WHERE bacterium = ?", (bacteriumName,))
	      print ''.ljust(60),bcolors.OKGREEN+'[ - CORRECTED - ]'+bcolors.ENDC
	      if args.log: logFile.write(bacteriumName+':\t REMOVED\n')


conn.commit()
conn.close() 
  
      
