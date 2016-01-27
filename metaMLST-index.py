import sys,os,subprocess,sqlite3,argparse,re,itertools
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


class bcolors:
	HEADER = '\033[95m'
	OKBLUE = '\033[94m'
	OKGREEN = '\033[92m'
	WARNING = '\033[93m'
	FAIL = '\033[91m'
	ENDC = '\033[0m'
	OKGREEN2 = '\033[42m\033[30m'

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

parser.add_argument("-d", "--database", help="database file")
parser.add_argument("-t", "--typings", help="typings in tab separated file (Build New Database)")
parser.add_argument("-s", "--sequences", help="Sequences (comma separated list of files) (Build New Database)")

parser.add_argument("-q","--dump_db", help="Dumps the database to a fasta file") 

parser.add_argument("-i","--buildindex", help="name of the output bowtie2 index") 
parser.add_argument("-b","--buildblast", help="name of the output blast index") 
parser.add_argument("-z","--buldblastfiler", help="filter on the blast index") 
args=parser.parse_args()


# Creates (or updates) database with -d -t -s options

if args.database and (args.typings or args.sequences):
	try:
		conn = sqlite3.connect(args.database)
	except IOError: 
		print "IOError: unable to access "+args.database+"!"

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
			print ('  ADDING SEQUENCES '+file).ljust(60)
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
					if re.match('^([a-zA-Z])*$',organism) and re.match('^([a-zA-Z0-9])*$',gene) and re.match('^([0-9])*$',allele):
						
						#if gene-allele couple is NOT present in database
						if len([row for row in cursor.execute("SELECT 1 FROM alleles WHERE bacterium = ? AND gene = ? and alleleVariant = ?",(organism,gene,allele))]) == 0:

							#Add gene-allele couple to the query-list for addition
							if gene not in [x for (x,k) in geneList]: geneList.append((gene,organism))
							alleleList.append((gene,organism,allele,str(sequence)))
							addCounter += 1
							
						else:
							print ('  Allele already present in DB: '+seq_record.id).ljust(60), bcolors.WARNING+'[ - Skip - ]'+bcolors.ENDC
							continue
					else:
						print ('  Invalid Sequence ID: '+seq_record.id).ljust(60), bcolors.FAIL+'[ - Skip - ]'+bcolors.ENDC
						continue
				else:
					print ('  Invalid Sequence ID: '+seq_record.id).ljust(60), bcolors.FAIL+'[ - Skip - ]'+bcolors.ENDC
					continue
					
			#todo aligned sequence
			#cursor.execute("SELECT 1 FROM genes WHERE bacterium = "); 
			print organism
			cursor.execute("INSERT OR IGNORE INTO organisms (organismkey) VALUES (?)",(organism,))
			cursor.executemany("INSERT OR IGNORE INTO genes (geneNAme, bacterium) VALUES (?,?)",geneList)
			cursor.executemany("INSERT INTO alleles (gene, bacterium,alleleVariant,sequence) VALUES (?,?,?,?)",alleleList)
			print (bcolors.OKGREEN+'[ - '+str(addCounter)+' PUSHED - ]'+bcolors.ENDC).rjust(15)
	
	if args.typings:
		for file in [seq.strip() for seq in args.typings.split(',')]:
			intest = 1
			profilesQuery=[]
			profilesLoaded=0
			fileOpen = open(file,'r')
			leng = int(len(fileOpen.readlines()))
			fileOpen.seek(0)
			profilesQuery = []
			
			for line in fileOpen:
				if line.startswith('@'): continue
				elif line == '': continue
				elif line.startswith('#'):
					organism = line.replace('#','').replace('_','').strip()
					
					print ('    DELETE old typings for : '+organism).ljust(60),
					sys.stdout.flush()
					cursor.execute("DELETE FROM profiles WHERE bacterium = ?",(organism,))
					print  bcolors.OKGREEN+'[ - DONE - ]'+bcolors.ENDC
					sys.stdout.flush()
					
					continue
					
				data = line.split()
				
				recID_Cache = dict((row['gene']+'_'+str(row['alleleVariant']),row['recID']) for row in cursor.execute("SELECT gene,alleleVariant,recID FROM alleles WHERE bacterium = ?",(organism,))) 
				
				
				if intest:
					print ('    READING MLST genes for : '+organism).ljust(60)
					sys.stdout.flush()
					intest = 0
					genes = data[1::]
					print ('    '+', '.join(data[1::])).ljust(60),bcolors.OKGREEN+'[ - DONE - ]'+bcolors.ENDC
					sys.stdout.flush()
				else:
					recIDs = []
					#print ('    PROFILE : '+data[0]).ljust(60)
					sys.stdout.flush()

					warningProfile=False #all Ok

					for key,variant in enumerate(data[1::]):
						
						if key < len(genes): 
							
							if (genes[key]+'_'+str(variant)) in recID_Cache: recIDs.append(recID_Cache[genes[key]+'_'+str(variant)])
							
							elif genes[key] in ['clonal_complex','species','mlst_clade']: continue
							else:
								print ('  Profile allele not found in DB : '+organism+'_'+genes[key]+'_'+variant).ljust(60), bcolors.FAIL+'[ - WARNING - ]'+bcolors.ENDC
								#sys.exit(1)  
								warningProfile=True
					 
					if not warningProfile: profilesQuery = itertools.chain(profilesQuery, [(organism,data[0],alleleRecID) for alleleRecID in recIDs])
					else: print ('  Profile Discarded: '+organism+'_'+genes[key]+'_'+variant).ljust(60), bcolors.FAIL+'[ - DONE - ]'+bcolors.ENDC

					print "\r    ANALYZING profile",organism,data[0],'|'+'-'*((profilesLoaded *10) / leng )+' '*(10-(profilesLoaded *10) / leng)+'|',str(round(float(profilesLoaded) / float(leng),4)*100)+'%',
					profilesLoaded +=1
					
					
			cursor.executemany("INSERT INTO profiles (bacterium, profileCode, alleleCode) VALUES (?,?,?)", profilesQuery)
					

			print ''.ljust(60),bcolors.OKGREEN+'[ - '+str(profilesLoaded)+' PUSHED - ]'+bcolors.ENDC
				
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
