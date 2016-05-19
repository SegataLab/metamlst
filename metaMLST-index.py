#!/usr/bin/env python

import sys,os,subprocess,sqlite3,argparse,re,itertools
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from metaMLST_functions import * 

def dump_db_to_fasta(path,dbPath):
	try:
		conn = sqlite3.connect(dbPath)
	except IOError: 
		print "IOError: unable to access "+dbPath+"!"

	conn.row_factory = sqlite3.Row
	cursor = conn.cursor() 
	
	print '  COLLECTING DATA'.ljust(60),
	sys.stdout.flush()
	
	strr=[SeqRecord(Seq(row['sequence'],IUPAC.unambiguous_dna),id=row['bacterium']+'_'+(row['gene']+'_'+str(row['alleleVariant'])),description='') for row in cursor.execute("SELECT bacterium,gene,alleleVariant,sequence FROM alleles WHERE sequence <> ''")]	
	SeqIO.write(strr,path,'fasta')
	print bcolors.OKGREEN+'[ - DONE - ]'+bcolors.ENDC

 
 
parser = argparse.ArgumentParser()

parser.add_argument("-d","--database", help="MetaMLST Database File (created with metaMLST-index", required=True)
parser.add_argument("-t", "--typings", help="typings in tab separated file (Build New Database)")
parser.add_argument("-s", "--sequences", help="Sequences (comma separated list of files)")

parser.add_argument("-q","--dump_db", help="Dumps the entire database to fasta file") 

parser.add_argument("-i","--buildindex", help="Build a Bowtie2 Index from the DB") 
parser.add_argument("-b","--buildblast", help="Build a BLAST Index from the DB") 
args=parser.parse_args()

try:
	conn = sqlite3.connect(args.database)
except IOError: 
	print "IOError: unable to access "+args.database+"!"
	sys.exit(0)


# Creates (or updates) database with -d -t -s options

if args.typings or args.sequences:

	conn.row_factory = sqlite3.Row
	cursor = conn.cursor() 
	
	cursor.execute("CREATE TABLE IF NOT EXISTS organisms (organismkey varchar(255), label VARCHAR(255), PRIMARY KEY(organismkey))")
	cursor.execute("CREATE TABLE IF NOT EXISTS genes (geneName varchar(255), bacterium VARCHAR(255), PRIMARY KEY(geneName,bacterium))")
	cursor.execute("CREATE TABLE IF NOT EXISTS alleles (recID INTEGER PRIMARY KEY AUTOINCREMENT,bacterium varchar(255), gene VARCHAR(255), sequence TEXT, alignedSequence TEXT, alleleVariant INT)")
	cursor.execute("CREATE TABLE IF NOT EXISTS profiles (recID INTEGER PRIMARY KEY AUTOINCREMENT, profileCode INTEGER, bacterium VARCHAR(255), alleleCode INTEGER)")
	
	# Alleles
	# Example:
	# >bacterium_recA_21
	# ATCTCTTGTGCT...
	# >recA22 
	
	
	if args.sequences:
		
		for file in [seq.strip() for seq in args.sequences.split(',')]:
			print (' ADDING SEQUENCES '+file+'...')
			alleleList = []
			geneList = []
			addCounter=0
			
			for seq_record in SeqIO.parse(file, "fasta"):
				splitLine = seq_record.id.split('_')
				
				if len(splitLine) == 3:
					allele = splitLine[2]
					gene = splitLine[1]
					organism = splitLine[0]
					sequence = seq_record.seq
					
					#IF line formatted in the correct way
					if re.match('^([a-zA-Z-])*$',organism) and re.match('^([a-zA-Z0-9|])*$',gene) and re.match('^([0-9])*$',allele):
						
						#if gene-allele couple is NOT present in database
						if len([row for row in cursor.execute("SELECT 1 FROM alleles WHERE bacterium = ? AND gene = ? and alleleVariant = ?",(organism,gene,allele))]) == 0:

							#Add gene-allele couple to the query-list for addition
							if gene not in [x for (x,k) in geneList]: geneList.append((gene,organism))
							alleleList.append((gene,organism,allele,str(sequence)))
							addCounter += 1
							
						else:
							print ('   Allele already present in DB: '+seq_record.id).ljust(50),(bcolors.FAIL+'[ - Skip - ]'+bcolors.ENDC).rjust(30) 
							continue
					else: 
						print ('   Invalid Sequence ID: '+seq_record.id).ljust(50),(bcolors.FAIL+'[ - Skip - ]'+bcolors.ENDC).rjust(30) 
						continue
				else:
					print ('   Invalid Sequence ID: '+seq_record.id).ljust(50),(bcolors.FAIL+'[ - Skip - ]'+bcolors.ENDC).rjust(30) 
					continue
					
			#todo aligned sequence
			#cursor.execute("SELECT 1 FROM genes WHERE bacterium = "); 
			
			cursor.executemany("INSERT OR IGNORE INTO genes (geneNAme, bacterium) VALUES (?,?)",geneList)
			cursor.executemany("INSERT INTO alleles (gene, bacterium,alleleVariant,sequence) VALUES (?,?,?,?)",alleleList)
			print (' ADDING SEQUENCES '+file).ljust(25), ('Added '+str(addCounter)+' seqs').rjust(25)+' '+(bcolors.OKGREEN+'[ - DONE - ]'+bcolors.ENDC).rjust(27) 
	
	if args.typings:
		for file in [seq.strip() for seq in args.typings.split(',')]:
			intest = 1
			profilesQuery=[]
			profilesLoaded=0
			fileOpen = open(file,'r')
			leng = int(len(fileOpen.readlines()))
			fileOpen.seek(0) 
			
			problematicList = {}
			for line in fileOpen:
				if line.startswith('@'): continue
				elif line == '': continue
				elif line.startswith('#'):
					if len(line.strip().split('|')) == 2:
						organism = line.strip().split('|')[0].replace('#','').replace('_','')
						organismextended = line.strip().split('|')[1]
					else: 
						organism = line.strip().split('|')[0].replace('#','').replace('_','')
						organismextended = line.strip().split('|')[0].replace('#','').replace('_','')
					
					cursor.execute("INSERT OR IGNORE INTO organisms (organismkey,label) VALUES (?,?)",(organism,organismextended))
					print (' DELETE old typings for '+organismextended).ljust(50),
					sys.stdout.flush()
					cursor.execute("DELETE FROM profiles WHERE bacterium = ?",(organism,))
					print  (bcolors.OKGREEN+'[ - DONE - ]'+bcolors.ENDC).rjust(30)
					sys.stdout.flush()
					
					continue
					     
				data = line.split()
				recID_Cache = dict((row['gene']+'_'+str(row['alleleVariant']),row['recID']) for row in cursor.execute("SELECT gene,alleleVariant,recID FROM alleles WHERE bacterium = ?",(organism,))) 
				problematic = False
				if intest:
					print (' READING MLST loci for '+organismextended).ljust(50)
					sys.stdout.flush()
					intest = 0
					genes = data[1::]
					print ('   '+', '.join([g for g in genes if g not in ['clonal_complex','species','mlst_clade']])).ljust(51)+(bcolors.OKGREEN+'[ - DONE - ]'+bcolors.ENDC).rjust(30) ### 
					sys.stdout.flush()
				else:
					recIDs = []
					sys.stdout.flush()
 
					
					for key,variant in enumerate(data[1::]):
						
						if key < len(genes): 
							
							if (genes[key]+'_'+str(variant)) in recID_Cache: recIDs.append(recID_Cache[genes[key]+'_'+str(variant)])
							
							elif genes[key] in ['clonal_complex','species','mlst_clade']: continue
							else:
								if str(data[0]) not in problematicList: problematicList[str(data[0])] = []
								problematicList[str(data[0])].append(organism+'_'+genes[key]+'_'+variant)
								problematic = True
				
					print "\r"+(" CHECKING PROFILES").ljust(24)+(' '+organism+' ST-'+data[0]).ljust(17) + ('[ - '+(str(int(float(profilesLoaded) / float(leng)*100))+'%').rjust(3)+' - ]').rjust(30),
					
					if not problematic:
						profilesLoaded +=1
						for element in recIDs:
							profilesQuery.append((organism,data[0],element))
							
						
			cursor.executemany("INSERT INTO profiles (bacterium, profileCode, alleleCode) VALUES (?,?,?)", profilesQuery)
			
			if len(problematicList) > 0:
				with open('metamlst_logfile.log','a') as logf:
					logf.write('The following STs for '+organism+' were skipped as one or more of alleles could not be found:\r\n')
					for key,element in problematicList.items():
						logf.write('ST-'+key+'\t'.join(element))
					logf.write(('-'*120)+'r\n')
				
			print '\r'+(' COMPLETED '+organism).ljust(26),('Added '+str(profilesLoaded)+' STs').rjust(25)+' '+(bcolors.OKGREEN+'[ - DONE - ]'+bcolors.ENDC).rjust(28) 
			
	conn.commit()
	conn.close() 
	# strr=[SeqRecord(Seq(row['sequence'],IUPAC.unambiguous_dna),id=row['bacterium']+'_'+(row['gene']+str(row['alleleVariant'])),description='') for row in cursor.execute("SELECT bacterium,gene,alleleVariant,sequence FROM alleles WHERE 1")]	
	# SeqIO.write(strr,'out.fa','fasta')

	# devnull = open('/dev/null', 'w')
	# child = subprocess.Popen("bowtie2-build out.fa out.index",shell=True, stdout=devnull)
	# child.wait()


if args.dump_db:
	dump_db_to_fasta(args.dump_db,args.database)
	
# Only builds BW2 index from DB (existing!)
if args.buildindex:
	
	dump_db_to_fasta('out.fa',args.database)

	print '  BUILDING INDEX'.ljust(60),
	sys.stdout.flush()
	
	devnull = open('/dev/null', 'w')
	child = subprocess.Popen("bowtie2-build out.fa "+args.buildindex,shell=True, stdout=devnull)
	child.wait()
	#os.remove('out.fa')
	print bcolors.OKGREEN+'[ - DONE - ]'+bcolors.ENDC
	conn.commit()
	conn.close()  

if args.buildblast:
	
	dump_db_to_fasta('out.fa',args.database)
	
	print '  BUILDING INDEX'.ljust(60),
	sys.stdout.flush()
	
	devnull = open('/dev/null', 'w')
	child = subprocess.Popen("makeblastdb -in out.fa -dbtype nucl -out "+args.buildblast,shell=True, stdout=devnull)
	child.wait()
	os.remove('out.fa')
	print bcolors.OKGREEN+'[ - DONE - ]'+bcolors.ENDC
	conn.commit()
	conn.close()  
