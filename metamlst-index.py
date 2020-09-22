#!/usr/bin/env python3

from __future__ import print_function

try:
	import sys,os,subprocess,sqlite3,argparse,re
	from metaMLST_functions import *
except ImportError as e:
	print("Error while importing python modules! Remember that this script requires: sys,os,subprocess,sqlite3,argparse,re")
	

	sys.exit(1)

try:
	from Bio import SeqIO
	from Bio.Seq import Seq
	from Bio.SeqRecord import SeqRecord
except ImportError as e:
	print("Error while importing Biopython. Please check Biopython is installed properly on your system!")
	sys.exit(1)


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		description='Builds and manages the MetaMLST SQLite Databases')

parser.add_argument("-t", "--typings", help="Typings in TAB separated file (Build New Database)")
parser.add_argument("-s", "--sequences", help="Sequences in FASTA format (comma separated list of files)")
parser.add_argument("-q","--dump_db", help="Dump the entire database to file in fasta format)") 
parser.add_argument("-i","--buildindex", help="Build a Bowtie2 Index from the DB") 
parser.add_argument("-b","--buildblast", help="Build a BLAST Index from the DB") 
parser.add_argument("-d",'--database', metavar="DB PATH", help="MetaMLST Database File (if unset, use the default database. If a file name is given, MetaMLST will create a new DB or update an existing one)")
parser.add_argument("--list", help="Lists all the MLST keys present in the database and exit", action="store_true") 
parser.add_argument("--filter", help="filters the db for a specific bacterium", default=None) 
parser.add_argument("--version", help="Prints version informations", action='store_true')

parser.add_argument("--bowtie2_threads", help="Number of Threads to use with bowtie2-build",default=4) 
parser.add_argument('--bowtie2_build', type=str, default='bowtie2-build',
        help="Full path to the bowtie2-build command to use, deafult assumes "
             "that 'bowtie2-build is present in the system path")


args=parser.parse_args()
if args.version:
	print_version()
	sys.exit(0)

try:
	#download the database if a non existing (but default-named) DB file is passed
	if not args.database:
		dbPath=check_install()
	else:
		dbPath=args.database

	metaMLSTDB = metaMLST_db(dbPath)
	conn = metaMLSTDB.conn
	cursor = metaMLSTDB.cursor

except IOError: 
	metamlst_print("Failed to connect to the database: please check your database file!",'FAIL',bcolors.FAIL)
	sys.exit(1)

if os.path.isfile(dbPath):
	try:
		cursor.execute("CREATE TABLE IF NOT EXISTS organisms (organismkey varchar(255), label VARCHAR(255), PRIMARY KEY(organismkey))")
		cursor.execute("CREATE TABLE IF NOT EXISTS genes (geneName varchar(255), bacterium VARCHAR(255), PRIMARY KEY(geneName,bacterium))")
		cursor.execute("CREATE TABLE IF NOT EXISTS alleles (recID INTEGER PRIMARY KEY AUTOINCREMENT,bacterium varchar(255), gene VARCHAR(255), sequence TEXT, alignedSequence TEXT, alleleVariant INT)")
		cursor.execute("CREATE TABLE IF NOT EXISTS profiles (recID INTEGER PRIMARY KEY AUTOINCREMENT, profileCode INTEGER, bacterium VARCHAR(255), alleleCode INTEGER)")
	
		print("Database "+dbPath+' contains:')
		cursor.execute("SELECT COUNT(*) as Mv FROM organisms WHERE 1")
		print ('\t'+str(cursor.fetchone()['Mv'])+' organisms')
		cursor.execute("SELECT COUNT(*) as Mv FROM genes WHERE 1")
		print ('\t'+str(cursor.fetchone()['Mv'])+' total loci')
		cursor.execute("SELECT COUNT(*) as Mv,SUM(LENGTH(sequence)) as Se FROM alleles WHERE 1")
		cont=cursor.fetchone()
		print ('\t'+str(cont['Mv'])+' total alleles (~'+str(round((cont['Se'] if cont['Se'] is not None else 0)/1000000.0,2))+' Megabases)')
		cursor.execute("SELECT COUNT(DISTINCT profileCode) as Mv FROM profiles WHERE 1")
		print ('\t'+str(cursor.fetchone()['Mv'])+' total profiles')

	except sqlite3.OperationalError as e:
		metamlst_print("Database Error: "+str(e),'FAIL',bcolors.FAIL)


if args.list:
	print ('-'*65)
	print ('MetaMLST Key'.ljust(30)+(' '*5)+'organism Full Name'.ljust(30))
	print ('-'*65)
	for key,label in db_getOrganisms(conn).items():
		print (key.ljust(30)+(' ')*5+label.ljust(30))
	sys.exit(0)

if args.typings or args.sequences:

	if args.sequences:
		
		for file in [seq.strip() for seq in args.sequences.split(',')]:
			metamlst_print('ADDING SEQUENCES','...',bcolors.HEADER)
			alleleList = []
			geneList = {}
			addCounter=0
			
			for seq_record in SeqIO.parse(file, "fasta"):
				splitLine = seq_record.id.split('_')
				
				if len(splitLine) == 3:
					allele = splitLine[2]
					gene = splitLine[1]
					organism = splitLine[0]
					sequence = seq_record.seq
					
					#IF line formatted in the correct way
					if re.match('^([a-zA-Z0-9-])*$',organism) and re.match('^([a-zA-Z0-9-])*$',gene) and re.match('^([0-9])*$',allele):
						
						#if gene-allele couple is NOT present in database
						if len([row for row in cursor.execute("SELECT 1 FROM alleles WHERE bacterium = ? AND gene = ? and alleleVariant = ?",(organism,gene,allele))]) == 0:

							#Add gene-allele couple to the query-list for addition
							
							#if gene geneList.append((gene,organism))
							if organism not in geneList: geneList[organism] = []
							if gene not in geneList[organism]: geneList[organism].append(gene)


							alleleList.append((gene,organism,allele,str(sequence)))
							addCounter += 1
							
						else:
							metamlst_print(' > Allele already present in DB: '+seq_record.id,'SKIP',bcolors.FAIL)
							continue
					else: 
						metamlst_print(' > Invalid Sequence ID: '+seq_record.id,'SKIP',bcolors.FAIL)
						continue
				else:
					metamlst_print(' > Malformed Sequence ID: '+seq_record.id,'SKIP',bcolors.FAIL)
					continue
			
			gListToAdd=[]		
			for org,genes in geneList.items():
				for gen in genes:
					gListToAdd.append((gen,org))
			
			cursor.executemany("INSERT OR IGNORE INTO genes (geneNAme, bacterium) VALUES (?,?)",gListToAdd)

			cursor.executemany("INSERT INTO alleles (gene, bacterium,alleleVariant,sequence) VALUES (?,?,?,?)",alleleList)
			metamlst_print('ADDING SEQUENCES '+file+' Added '+str(addCounter)+' seqs','DONE',bcolors.OKGREEN)
	
	if args.typings:
		for file in [seq.strip() for seq in args.typings.split(',')]:
			intest = 1
			profilesQuery=[]
			profilesLoaded=0
			fileOpen = open(file,'r')
			leng = int(len(fileOpen.readlines()))-2
			fileOpen.seek(0) 
			
			problematicList = {}
			for line in fileOpen:
				if line.startswith('@'): continue
				elif line == '': continue
				elif line.startswith('#'):
					organism = line.strip().split('|')[0].replace('#','').replace('_','')
					organismLabel = line.strip().split('|')[1] if len(line.strip().split('|')) == 2 else organism
					
					cursor.execute("INSERT OR IGNORE INTO organisms (organismkey,label) VALUES (?,?)",(organism,organismLabel))
					cursor.execute("DELETE FROM profiles WHERE bacterium = ?",(organism,))
					metamlst_print('DELETED all profiles for '+organismLabel,'DONE',bcolors.OKGREEN)
					 
					continue
					     
				data = line.split()
				recID_Cache = dict((row['gene']+'_'+str(row['alleleVariant']),row['recID']) for row in cursor.execute("SELECT gene,alleleVariant,recID FROM alleles WHERE bacterium = ?",(organism,))) 
				problematic = False
				if intest:
					
					metamlst_print('READING MLST loci for '+organismLabel,'....',bcolors.OKGREEN)
					sys.stdout.flush()
					intest = 0
					genes = data[1::]
					
					metamlst_print(' > '+', '.join([g for g in genes if g not in MLST_KEYWORDS]),'DONE',bcolors.OKGREEN)
					
					sys.stdout.flush()
				else:
					recIDs = []
					sys.stdout.flush()
 
					
					for key,variant in enumerate(data[1::]):
						
						if key < len(genes): 
							
							if (genes[key]+'_'+str(variant)) in recID_Cache: recIDs.append(recID_Cache[genes[key]+'_'+str(variant)])
							
							elif genes[key] in ['clonal_complex','clonal-complex','species','mlst_clade']: continue
							else:
								if str(data[0]) not in problematicList: problematicList[str(data[0])] = []
								problematicList[str(data[0])].append(organism+'_'+genes[key]+'_'+variant)
								problematic = True
				
					#print "\r"+(" CHECKING PROFILES").ljust(24)+(' '+organism+' ST-'+data[0]).ljust(17) + ('[ - '+(str(int(float(profilesLoaded) / float(leng)*100))+'%').rjust(3)+' - ]').rjust(30),
					percentCompleted = float(profilesLoaded) / float(leng)*100.0

					metamlst_print('CHECKING PROFILE ['+organism+ ' ' + data[0] + ']',str(int(percentCompleted))+'%',bcolors.OKGREEN,reline=True)

					#print (problematicList)
					if not problematic:
						profilesLoaded +=1
						for element in recIDs:
							profilesQuery.append((organism,data[0],element))
			if profilesLoaded>0:
				cursor.execute("INSERT OR IGNORE INTO organisms (organismkey,label) VALUES (?,?)",(organism,organismLabel))

			#last update to percentage
			percentCompleted = float(profilesLoaded) / float(leng)*100.0

			metamlst_print(str(profilesLoaded)+'/'+str(leng)+' PROFILES LOADED',str(int(percentCompleted))+'%',bcolors.OKGREEN,reline=True,newLine=True)
			
			cursor.executemany("INSERT INTO profiles (bacterium, profileCode, alleleCode) VALUES (?,?,?)", profilesQuery)
			
			if len(problematicList) > 0:
				with open('metamlst_logfile.log','a') as logf:
					logf.write('The following STs for '+organism+' were skipped as one or more of the alleles comprising the profile could not be found in your DB:\r\n')
					for key,element in problematicList.items():
						logf.write('ST-'+' '+key+'\t'.join(element)+' was missing \r\n')
					logf.write(('-'*120)+'r\n')
				
			#print '\r'+(' COMPLETED '+organism).ljust(26),('Added '+str(profilesLoaded)+' STs').rjust(25)+' '+(bcolors.OKGREEN+'[ - DONE - ]'+bcolors.ENDC).rjust(28) 
			metamlst_print('COMPLETED '+organism,'DONE',bcolors.OKGREEN)

if args.dump_db:
	dump_db_to_fasta(conn,args.dump_db,args.filter)
	
if args.buildindex:

	if os.path.isfile(args.buildindex+'.1.bt2'):
		metamlst_print("File "+args.buildindex+'.1.bt2 is already present. It seems the index is already there.','SKIP',bcolors.FAIL)
	else:
		dump_db_to_fasta(conn,'out.fa',args.filter)
		metamlst_print('BUILDING INDEX','...',bcolors.HEADER)
		sys.stdout.flush()
		

		bt2_cmd = [args.bowtie2_build, '--quiet']

		if args.bowtie2_threads > 1:
			bt2_build_output = subprocess.check_output([args.bowtie2_build  , '--usage'], stderr=subprocess.STDOUT)

			if 'threads' in str(bt2_build_output):
				bt2_cmd += ['--threads', str(args.bowtie2_threads)]

		bt2_cmd += ['-f', 'out.fa', args.buildindex]

		try:
			subprocess.check_call(bt2_cmd)
			os.remove('out.fa')
			metamlst_print('BUILDING INDEX','DONE',bcolors.OKGREEN,reline=True,newLine=True) 
		except Exception as e:
			metamlst_print('Fatal error running bowtie2-index. Error message: '+e,'!',bcolors.FAIL) 
			sys.exit(1)

	
if args.buildblast:
	
	dump_db_to_fasta(conn,'out.fa',args.filter)
	
	metamlst_print('BUILDING INDEX','...',bcolors.HEADER)
	sys.stdout.flush()
	
	devnull = open('/dev/null', 'w')
	child = subprocess.Popen("makeblastdb -in out.fa -dbtype nucl -out "+args.buildblast,shell=True, stdout=devnull)
	child.wait()
	metamlst_print('BUILDING INDEX','DONE',bcolors.OKGREEN,reline=True,newLine=True)

conn.commit()
metaMLSTDB.closeConnection()  
