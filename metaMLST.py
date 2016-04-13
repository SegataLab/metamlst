#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import json, math,sqlite3,itertools,subprocess,time,gc,argparse, re, os,sys
from metaMLST_functions import *



parser = argparse.ArgumentParser()
parser.add_argument("file", help="BAM file containing the sequences")
parser.add_argument("--filter", help="focus on specific set of organisms only (METAMLST-KEY, space separated)",  nargs='+')

parser.add_argument("-p","--penalty", help="penalty for not found allele", default=100, type=int)
parser.add_argument("--minscore", help="Minimum BowTie2 score (absolute val!)", default=80, type=int)

#parser.add_argument("--silent", help="Silent Mode", action='store_true')
parser.add_argument("--debug", help="Debug Mode", action='store_true', default=False) 
parser.add_argument("--out_folder", help="Output Folder", default='./out')

parser.add_argument("--mincoverage", help="minimum coverage to produce FASTQ assembler files", default=5, type=int)
parser.add_argument("--min_read_len", help="minimum length for bowtie2 reads", default=50, type=int)
parser.add_argument("--max_xM", help="maximum XM", default=5, type=int)
parser.add_argument("--min_accuracy", help="minimum percentage of allele-reconstruction for each organism to pass", default=0.90, type=float)

parser.add_argument("--present_loci", help="Percentage of MLST loci to type an organism, default: 100", default=100, type=int)
parser.add_argument("--log", help="generate logfiles", action="store_true") 
parser.add_argument("-d","--database", help="MetaMLST Database File (created with metaMLST-index", required=True)
parser.add_argument("-a","--all_sequences", help="output also known sequences", action="store_true")

args=parser.parse_args()

#PREPARE
MetaMLSTDBconn = sqlite3.connect(args.database)
MetaMLSTDBconn.row_factory = sqlite3.Row
c = MetaMLSTDBconn.cursor()

fil = open(args.file,'r')
cel = {}
geneMaxiMin = {}
sequenceBank = {}
ignoredReads = 0
 

fileName = (args.file.split('/'))[-1].split('.')[0]

if not os.path.isdir(args.out_folder): os.mkdir(args.out_folder)
workUnit = args.out_folder	
 
child = subprocess.Popen("samtools view -h - ",shell=True, stdout=subprocess.PIPE, stdin = fil)
 
#samFileContent = fil.readlines()

for line in child.stdout:
	if(line[0] != '@'): 
		read = line.split('\t')
		readCode = read[0]
		# species = re.findall('[a-zA-Z]*_',read[2])[0].replace('_','')
		# gene = re.findall('_[a-zA-Z_]*',read[2])[0].replace('_','')
		# allele = re.findall('[0-9]*$',read[2])[0]
		
		# (Added)
		species,gene,allele = read[2].split('_')
		# (Removed)
		# species = read[2].split('_')[0]
		# gene = read[2].split('_')[1]
		# allele = read[2].split('_')[2]
		
		score = (int)((read[11]).split(':')[2])
		xM = (int)((read[14]).split(':')[2])
		sequence = read[9]
		quality = read[10]
		#quality = [(ord(x)-ord("!")) for x in read[10]] ???
		 
		if (args.filter and species in args.filter) or not args.filter:
			if score >= args.minscore and len(sequence) >= args.min_read_len and xM <= args.max_xM:
				if species not in cel: cel[species] = {} #new species
			
				if gene not in cel[species]: #found a new gene
					cel[species][gene] = {}
					if args.log: geneMaxiMin[species+gene] = [0,float("inf")]
					#Removed: (De-comment to have sequences and quality for FASTQ generation)
					sequenceBank[species+'_'+gene] = {}
					# sequenceBank[species+'_'+gene] = []
					
				if allele not in cel[species][gene]: #found a new allele
					cel[species][gene][allele] = []
				
				cel[species][gene][allele].append(score)
				#sequenceBank[species+'_'+gene].append(SeqRecord(Seq(sequence, IUPAC.unambiguous_dna), id = readCode, description = ''))
				#sequenceBank[species+'_'+gene].append((readCode,sequence,quality))
				
				# Removed: (De-comment to have seuqences and quality for FASTQ generation)
					# sequenceBank[species+'_'+gene][readCode]=(sequence,quality)
				sequenceBank[species+'_'+gene][readCode]=len(sequence)
				# sequenceBank[species+'_'+gene][readCode] 
				
				#print readCode
			else: ignoredReads = ignoredReads+1

#COMPILES THE DATA STRUCTURE 			
for speciesKey,species in cel.items():
	for geneKey, geneInfo in species.items(): #geni
	
		maxLen = max([len(x) for x in geneInfo.values()])

		for geneInfoKey,geneValues in geneInfo.items(): #alleli passata 2
			geneLen = len(geneValues)
			
			localScore = sum(item for item in geneValues)
			
			
			
			if geneLen != maxLen:
				localScore = localScore - (maxLen-geneLen)*args.penalty
			averageScore = float(localScore)/float(geneLen)
			cel[speciesKey][geneKey][geneInfoKey] = (localScore,maxLen,round(averageScore,1)) 
		
			if args.log:
				if cel[speciesKey][geneKey][geneInfoKey][0] > geneMaxiMin[speciesKey+geneKey][0]:
					geneMaxiMin[speciesKey+geneKey][0] = cel[speciesKey][geneKey][geneInfoKey][0]
					
				if cel[speciesKey][geneKey][geneInfoKey][0] < geneMaxiMin[speciesKey+geneKey][1]:
					geneMaxiMin[speciesKey+geneKey][1] = cel[speciesKey][geneKey][geneInfoKey][0]


		

# cel ['haemophilus']['haemophilus_gene'][allele] = (a, b, c) --> a: summed score, b = maxLen, c = a/b

###################################### OUT FILE 1 (LOG)
if args.log:
	dfil = open(workUnit+'/'+fileName+'_'+str(int(time.time()))+'.out','w')
	dfil.write("SAMPLE:\t\t\t\t\t"+args.file+'\r\n')
	dfil.write("PENALTY:\t\t\t\t"+repr(args.penalty)+'\r\n')
	dfil.write("MIN-THRESHOLD SCORE:\t"+repr(args.minscore)+'\r\n')
	dfil.write("IGNORED:\t\t\t\t"+repr(ignoredReads)+' BAM READS\r\n\r\n------------------------------  RESULTS ------------------------------\r\n')

	for speciesKey,species in cel.items():
		dfil.write(speciesKey)
		dfil.write('\r\n{')
		for geneKey, geneInfo in species.items(): #geni
			dfil.write('\t'+geneKey+'\t'+repr((geneMaxiMin[speciesKey+geneKey])[::-1])+'\r\n')
			dfil.write('\t{\r\n')
			for geneInfoKey,(score,maxLen,average) in sorted(geneInfo.items(), key= lambda x: x[1]): #alleli ordinati
				dfil.write('\t\t'+geneKey + geneInfoKey+'\t\t'+repr(score)+'\t\t'+repr(maxLen)+'\t\t'+str(average)+'\r\n')
			dfil.write('\t}\r\n')
		dfil.write('}\r\n')
	dfil.close()	


 
print '\r\n'+bcolors.OKBLUE+('  '+fileName+'  ').center(75,'-')+bcolors.ENDC

for speciesKey,species in cel.items():
	
	#tVar is the dictionary of locus detected in the organism
	tVar = dict([(row['geneName'],0) for row in  c.execute("SELECT geneName FROM genes WHERE bacterium = ?",(speciesKey,))])
	
	#GENE PRESENCE 
	
	if len(tVar) < len(species.keys()):
		print 'Database is broken for' +speciesKey+bcolors.FAIL+'[ - EXITING - ]'.rjust(75,' ')+bcolors.ENDC
		sys.exit(0)
	
	for sk in species.keys():
		tVar[sk] = 1
	vals = sum([t for t in tVar.values()]) 
	
	
	print ((bcolors.OKGREEN if (int((float(vals)/float(len(tVar)))*100) >= args.present_loci) else bcolors.FAIL)+' '+speciesKey.ljust(18,' ')+bcolors.ENDC)+' Detected Loci: '+', '.join([bcolors.OKGREEN + sk + bcolors.ENDC for (sk,v) in sorted(tVar.items(), key = lambda x:x[0]) if v == 1])
	print (' '*20)+'Missing Loci : '+', '.join([bcolors.FAIL + sk + bcolors.ENDC for (sk,v) in sorted(tVar.items(), key = lambda x:x[0]) if v == 0])
	
		
	if int((float(vals)/float(len(tVar)))*100) >= args.present_loci:
	
		print " \r\n >> Closest allele identification <<"
		
		print "\r\n  "+"Locus".ljust(7)+"Avg. Coverage".rjust(15)+" Closest Allele(s)".ljust(36)+"Score".rjust(7)+"Hits".rjust(5)
		sys.stdout.flush()
		
		for geneKey, geneInfo in sorted(species.items(),key=lambda x:x[0]): #loci
		
			minValue = max([avg for (val,leng,avg) in geneInfo.values()])
			
			aElements = {}
			tmp=[]
			for k,(val,leng,avg) in geneInfo.items():
				if avg == minValue:
					aElements[k]=(val,leng,avg)
					tmp.append(k)
			tmp = ",".join(sorted(tmp))
			
			sequenceKey = speciesKey+'_'+geneKey
			c.execute("SELECT LENGTH(sequence) as L FROM alleles WHERE bacterium = ? AND gene = ? ORDER BY L DESC LIMIT 1", (speciesKey,geneKey))
			genL = c.fetchone()['L']
			#Removed: de-comment to have sequences and quality for FASTQ generation
				# coverage = sum([len(x) for (x,q) in sequenceBank[sequenceKey].values()])
			# print [x for (x,q) in sequenceBank[sequenceKey].values()]
			coverage = sum([x for x in sequenceBank[sequenceKey].values()])
			
			# if coverage >= args.mincoverage * genL: 
				# color = "\033[92m"
				# Removed: de-comment to have sequences and quality for FASTQ generation
					# fqfil = open(workUnit+'/'+sequenceKey+'.fasta','a')
					# for sequenceSpec,(sequence,quality) in sequenceBank[sequenceKey].items():
						# fqfil.write('>'+sequenceSpec+'\r\n')  
						# fqfil.write(sequence+'\r\n')   
					# fqfil.close()
			# else:  color = "\033[93m"
			
			print "  "+bcolors.WARNING+geneKey.ljust(7)+bcolors.ENDC+str(round(float(coverage)/float(genL),2)).center(15)+bcolors.ENDC+bcolors.OKBLUE,tmp.ljust(36)+bcolors.ENDC+bcolors.HEADER+str(minValue).rjust(6)+str(aElements.itervalues().next()[1]).rjust(5)+bcolors.ENDC
		print ""
		sys.stdout.flush()
		# matchingProfiles = []
		# for row in c.execute("SELECT profileCode, COUNT(*) as T FROM profiles WHERE alleleCode IN ("+','.join(profileTrack)+") GROUP BY profileCode HAVING T = (SELECT COUNT(*)  FROM profiles WHERE alleleCode IN ("+','.join(profileTrack)+") GROUP BY profileCode ORDER BY COUNT(*) DESC LIMIT 1) ORDER BY T DESC"):
			# matchScore = str(round(float(row['T']) / float(len(tVar)),4)*100)+' %'
			# print ("  MLST PROFILE "+str(row['profileCode'])).ljust(23)+str("MATCH: "+matchScore).rjust(10)
			# matchingProfiles.append(row['profileCode'])
			#v = [riw['sequence'] for riw in ]
		
		
		sampleSequence = '' #for organism
		
		
		# for geneKey, geneInfo in sorted(species.items(), key= lambda x: x[0]):
			# alleles = [geneKey+'_'+str(k) for k,(val,leng,avg) in geneInfo.items() if avg == min([avg1 for (val1,leng1,avg1) in geneInfo.values()])]
			# for allele in alleles:
		
		print " >> Consensous Sequence Build <<"
			
		l = [sorted([(speciesKey+'_'+g1+'_'+k,db_getUnalSequence(MetaMLSTDBconn,speciesKey,g1,k)) for k,(val,leng,avg) in g2.items() if avg == max([avg1 for (val1,leng1,avg1) in g2.values()])])[0] for g1,g2 in species.items()] 
		consenSeq = buildConsensus(args.file, dict(l),args.minscore,args.max_xM,args.debug)
		

		newProfile = 0
		finWrite = 1
		#print "      "+"Gene".ljust(6)+'Ref.'.ljust(7)+"Length".rjust(10)+"Ns".rjust(10)+"SNPs".rjust(10)
		print "\r\n  "+"Gene".ljust(6)+"Ref.".ljust(7)+"Length".rjust(7)+"Ns".rjust(7)+"SNPs".rjust(7)+"Breadth of Coverage".rjust(21)+"Notes".rjust(10)
		for l in consenSeq:
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
				if not args.all_sequences: l.seq = ''
			
			print "  "+(l.id.split('_')[1]).ljust(6)+(l.id.split('_')[2]).ljust(7)+leng.rjust(7)+holes.rjust(7)+str(snps).rjust(7)+leng_ns.rjust(21)+newAllele.rjust(10)
			
		
		print ''
		 
		
		if finWrite:
			print '  Reconstruction Successful > Write'
			profil = open(workUnit+'/'+fileName+'.nfo','a')	 
			profil.write(speciesKey+'\t'+fileName+'\t'+"\t".join([recd.id+"::"+str(recd.seq)+'::'+str(round(1-float(recd.description.split('_')[0].split('::')[1])/float(recd.seqLen),4)*100)+'::'+str(round(float(recd.description.split('_')[1].split('::')[1])/float(recd.seqLen),4)*100) for recd in consenSeq])+'\r\n')
			
			profil.close()
		else: print '  Reconstruction Failed  ('+str(round(args.min_accuracy*100,2))+'% min.accuracy) > \033[91mSkip\033[0m'
		print ''
		
	del cel[speciesKey]
	gc.collect()
			
 
MetaMLSTDBconn.close() 

if len(cel): print '\033[92m'+'[ - Completed - ]'.rjust(75,' ')+'\033[0m'

