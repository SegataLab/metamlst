#!/usr/bin/env python2.7

import sys,os,subprocess,argparse,os

from metaMLST_functions import *

try:
    from Bio import SeqIO
except ImportError:
    sys.stderr.write("Error! Biopython SeqIO library not detected!\n")
    sys.exit(1)

try:
    from Bio.Seq import Seq
except ImportError:
    sys.stderr.write("Error! Biopython Seq library not detected!\n")
    sys.exit(1)

try: from Bio.SeqRecord import SeqRecord
except ImportError:
    sys.stderr.write("Error! Biopython SeqRecord library not detected!\n")
    sys.exit(1)
    
try: from Bio.Alphabet import IUPAC 
except ImportError:
    sys.stderr.write("Error! Biopython Alphabet library not detected!\n")
    sys.exit(1)
    
try: from cStringIO import StringIO
except ImportError:
    sys.stderr.write("Error! Biopython cStringIO library not detected!\n")
    sys.exit(1)


parser = argparse.ArgumentParser('Performs MLST analysis on contigs or genomes')
parser.add_argument("files", help="Input .fasta files (can be a folder)",default="")
parser.add_argument("-d",'--database', metavar="DB PATH", help="Specify a different MetaMLST-Database. If unset, use the default Database. You can create a custom DB with metaMLST-index.py)")
parser.add_argument("-w","--work", help="Output files will be placed in this folder. By default the current folder is used (./) ", default='.')
parser.add_argument("--quiet", help="No output on stdin", action="store_true")
parser.add_argument("--min_pident", help="Minimum percentage of identity to the reference for each each BLAST to be consideretd (default: 90)", default=90.0, type=float)
parser.add_argument("--min_length", help="Minimum percentage of gene-coverage for each BLAST alignment to be considered (default: 90)", default=90.0, type=float)
parser.add_argument("--blastdb_prefix", help="Overrides the creation of a BLASTDB, and use a custom one")
parser.add_argument("--version", help="Prints version informations", action='store_true')
parser.add_argument("profile", help="MLST key (e.g. ecoli). To see all the available profiles use '?' as profile", default="")
args=parser.parse_args()

if args.version: print_version()
if args.files == '':
	parser.print_help()
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

if args.profile=='?':
	print 'Organism Name'.ljust(30)+(' '*5)+'MetaMLST key'.ljust(30)
	print '-'*65
	for key,label in metaMLSTDB.getOrganisms().items():
		print key.ljust(30)+(' ')*5+label.ljust(30)
	sys.exit(0)

masterLog = open(args.work+'/data_'+args.profile+'.txt','w')

if not args.blastdb_prefix: #if no index, create
		
	seqStream = StringIO()
	SeqIO.write(metaMLSTDB.getAlleles(args.profile),seqStream,'fasta') 
		
	devnull = open('/dev/null', 'w')
	child = subprocess.Popen("makeblastdb -in - -dbtype nucl -title abc -out mlstdbx_"+args.profile,shell=True, stdout=devnull, stdin=subprocess.PIPE)
	child.communicate(seqStream.getvalue())
	child.wait() 

profileKeys=metaMLSTDB.getGeneNames(args.profile)
masterLog.write('SAMPLE\tBACTERIUM\tST\tST_ACCURACY\t'+'\t'.join([k+'\t'+k+'_perc_iden\t'+k+'_len_of_gene\t'+k+'_len_aligned' for k in sorted(profileKeys)])+'\r\n')

if not args.quiet:
	print(bcolors.OKGREEN+"Long/Exact Match"+bcolors.ENDC+'\t No SNPs: perfect match')
	print(bcolors.WARNING+"Short/Exact Match "+bcolors.ENDC+'\t No SNPs, part of locus is covered  (report closest)')
	print(bcolors.OKBLUE+"Long/Partial Match"+bcolors.ENDC+'\t Some SNPs, whole locus is covered (report closest)')
	print(bcolors.RED+"Short/Partial Match"+bcolors.ENDC+'\t Some SNPs, part of locus is covered (report closest)\n\n')

prefix = ''
if os.path.isdir(args.files):
	prefix=args.files+'/'
	subFiles = os.listdir(args.files)
else: subFiles = args.files.split(',')
	
for file in subFiles:
	
	if file.split('.')[-1] not in ['fa','fss','ffn','fasta','fna','faa']: continue 
	scor = {}
	blasted= []
	dbb = 'mlstdbx_'+args.profile if not args.blastdb_prefix else args.blastdb_prefix
	
	child = subprocess.Popen('blastn -query '+prefix+file+' -max_target_seqs 100000 -outfmt "6 qseqid sseqid qlen slen length pident qseq sseq sstart send score" -db '+dbb,shell=True, stdout=subprocess.PIPE, stdin = subprocess.PIPE)
	out = StringIO(child.communicate()[0]) 

	for line in out.getvalue().split('\n'):  
		lat = line.split('\t')	
		if len(lat) < 4: continue
		(organism,gene,allele) = lat[1].split('_')
		
		if gene not in scor: scor[gene] = []
		entry = {'label':lat[1],'slen':lat[3],'leng':lat[4],'pident':lat[5],'score':int(lat[10]),'seq':lat[6],'start':lat[8],'end':lat[9]}
		if float(entry['pident']) >= args.min_pident and ((float(entry['leng'])/float(entry['slen']))*100) >= args.min_length: scor[gene].append(entry)
		
	for gen,listval in scor.items():
	
		if (len(listval) == 0 ): continue
		maxScore = max([element['score'] for element in listval])
		aligns = sorted([element for element in listval if element['score'] == maxScore],key=lambda x: x['label'])[0]		
		# aligns = sorted([element for element in listval if element['score'] == maxScore],key=lambda x : float(x['leng']) / float(x['slen']) * float(x['pident']))[0]
		#else: aligns = [element for element in sorted(listval,reverse=True,key=lambda x : float(x['leng']) / float(x['slen']) * float(x['pident']))][0]
		 
		blasted.append( aligns )
	

	allelic = dict((kel,{}) for kel in profileKeys)
	profilic = []
	pstrings = {}
	pseqs = {}
	
	#qseqid sseqid qlen slen length pident qseq sseq
	
	for cdict in blasted: #for each entry (gene)
		target = cdict['label']
		perc = float(cdict['pident'])
		slen = int(cdict['slen']) #length of the gene (subject = DB)
		leng = int(cdict['leng']) #length of the match (genomic read)
		qSeq = Seq(cdict['seq'],IUPAC.unambiguous_dna)
		# qSeq =   #length of the match (genomic read)
		qstart = int(cdict['start'])
		qend = int(cdict['end'])
		
		(organism,gene,allele) = target.split('_')
		  
		if qstart > qend:#reverse
			qSeq = qSeq.reverse_complement()
			qstart,qend = qend,qstart 
		   
		dashSequence = SeqRecord(Seq('-'*(qstart-1)+str(qSeq)+'-'*(slen-qend),IUPAC.unambiguous_dna),id=target+'_'+str(perc)+'_'+str(leng)+'/'+str(slen),description='')		
		
		#Correct Len, Perfect match			GREEN
		#Correct Len, Non-perfect match		YELLOW 
		#Wrong Len, Perfect match			CYAN
		#Wrong Len, Non-perfect match		RED
		if perc == 100.0 and slen == leng: color = (bcolors.OKGREEN,'','','')
		elif perc == 100.0 and slen != leng: color = (bcolors.WARNING,'',bcolors.WARNING,'*')
		elif perc != 100.0 and slen == leng: color = (bcolors.OKBLUE,bcolors.OKBLUE,'','*')
		elif perc != 100.0 and slen != leng: color = (bcolors.RED,bcolors.RED,bcolors.RED,'*')
		 
		allelic[gene]['allele'] = str(allele)
		allelic[gene]['perc'] = str(perc)
		allelic[gene]['len'] = str(leng)+'/'+str(slen)
		allelic[gene]['leng'] = str(leng)
		allelic[gene]['slen'] = str(slen)
		allelic[gene]['color'] = color  
		allelic[gene]['target'] = target  
		
		allelic[gene]['sequence'] = dashSequence  
		
		if perc == 100.0:  profilic.append(target)
		
	
	
	profileID = '--'
	profileScore = '--'
	if len([1 for v in allelic.values() if v == {}]) == 0:
		
		tryDefine = metaMLSTDB.defineProfile(profilic)

		if tryDefine:
			profileID = str(tryDefine[0][0])
			profileScore = str(tryDefine[0][1])
			
			#metaMLST-like-file phase
			af = open(args.work+'/'+os.path.basename(file).replace('.fna','')+'.nfo','a') 
			af.write( args.profile+'\t'+os.path.basename(file)+'\t'+'\t'.join([str(allelicElement['target'])+'::'+(str(allelicElement['sequence'].seq) if ( float(allelicElement['leng'])/float(allelicElement['slen']) * float(allelicElement['perc']) != 100.0 ) else '')+'::'+'100.0::0.0' for allKey,allelicElement in sorted(allelic.items()) if allelicElement != {}])+'\r\n')
			af.close()	
			#print os.path.basename(file)+'\t'+profileID+'\t'+profileScore
			
	

	#Output Phase
	
	synthFile = StringIO()
	SeqIO.write([allelicElement['sequence'] for allKey,allelicElement in sorted(allelic.items()) if allelicElement != {}],synthFile,'fasta')
	of = open(args.work+'/report_'+os.path.basename(file)[:15]+'.txt','w')
	of.write('\n\n#TABLE OF RESULTS: '+os.path.basename(file)+'\n\n')
	of.write( '#'+''.rjust(18)+''.join([k.center(11) for k in sorted(profileKeys)])+'ST'.center(11)+'\n' )
	of.write( '#'+'Allelic Profile'.rjust(18)+''.join([ (allelic[k]['allele']+allelic[k]['color'][3]).center(11) if allelic[k] != {} else '-'.center(11) for k in sorted(allelic.keys())])+(profileID+' ('+str(profileScore)+'%)').center(14)+'\n')
	of.write( '#'+'Perc. Ident.'.rjust(18)+''.join([allelic[k]['perc'].center(11)   if allelic[k] != {} else '-'.center(11) for k in sorted(allelic.keys())])+'|\n' )
	of.write( '#'+'Length.'.rjust(18)+''.join([allelic[k]['len'].center(11) if allelic[k] != {} else '-'.center(11) for k in sorted(allelic.keys())])+'|\n' )
	
	of.write('\n\n#SEQUENCES\n\n'+synthFile.getvalue())
	of.close()
	
	if not args.quiet:
		print ('FILE'.ljust(15)+'|'.join([k.center(7) for k in sorted(profileKeys)])+'|'+'ST'.center(5))
		print ('-'*80)
		print (os.path.basename(file)[:14].ljust(15)+'|'.join([ (allelic[k]['color'][0]+allelic[k]['allele'].center(7)+bcolors.ENDC) if allelic[k] != {} else '-'.center(7) for k in sorted(allelic.keys())])+'|'+(profileID+' ('+str(profileScore)+'%)').center(7))
		print ('-'*80)
		print ('Perc. Ident.  '.rjust(15)+'|'.join([ (allelic[k]['color'][1]+allelic[k]['perc'].center(7)+bcolors.ENDC)   if allelic[k] != {} else '-'.center(7) for k in sorted(allelic.keys())])+'|')
		print ('Length.  '.rjust(15)+'|'.join([(allelic[k]['color'][2]+allelic[k]['len'].center(7)+bcolors.ENDC) if allelic[k] != {} else '-'.center(7) for k in sorted(allelic.keys())])+'|')
		print ('')
	
	masterLog.write(os.path.basename(file)+'\t'+args.profile+'\t'+profileID+'\t'+profileScore+'\t'+'\t'.join( [(allelic[k]['allele']+'\t'+allelic[k]['perc']+'\t'+allelic[k]['leng']+'\t'+allelic[k]['slen']) if allelic[k] != {} else '-\t-\t-\t-' for k in sorted(allelic.keys())] )+'\r\n' )
	
	
masterLog.close()

if not args.blastdb_prefix:
	os.remove('mlstdbx_'+args.profile+'.nin')
	os.remove('mlstdbx_'+args.profile+'.nsq')
	os.remove('mlstdbx_'+args.profile+'.nhr')
