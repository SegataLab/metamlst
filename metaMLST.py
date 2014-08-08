import argparse, re, os,sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from StringIO import StringIO
import json, math,sqlite3,itertools,subprocess,time,gc
 
def alleleInProfile(bacterium,gene,allele,profile):
	e = conn.cursor()
	e.execute("SELECT 1 FROM profiles,alleles WHERE alleleCode = alleles.recID AND profileCode = ? AND alleleVariant = ? AND gene = ? AND profiles.bacterium = ?",(profile,allele,gene,bacterium))
	return (len(e.fetchall()) > 0)
	
def sequenceExists(bacterium,sequence):
	
	e = conn.cursor()
	e.execute("SELECT 1 FROM alleles WHERE sequence = ? AND bacterium = ?",(str(sequence),bacterium))
	return (len(e.fetchall()) > 0)

def sequenceFind(bacterium,sequence):
	
	e = conn.cursor()
	e.execute("SELECT gene,alleleVariant FROM alleles WHERE sequence = ? AND bacterium = ?",(str(sequence),bacterium))
	res = e.fetchone()
	
	if res: return res['gene']
	else: return 0

def alleleInManyProfile(bacterium,gene,allele,profileList):
	
	e = conn.cursor()
	profileString = ','.join([str(x) for x in profileList])
	e.execute("SELECT 1 FROM profiles,alleles WHERE alleleCode = alleles.recID AND profileCode IN ("+profileString+")  AND alleleVariant = ? AND gene = ? AND profiles.bacterium = ?",(allele,gene,bacterium))
	return (len(e.fetchall()) > 0)
	
def db_getSequence(bacterium,gene,allele):
	e = conn.cursor()
	e.execute("SELECT alignedSequence FROM alleles WHERE bacterium = ? AND gene = ? AND alleleVariant = ?",(bacterium,gene,allele))
	return (e.fetchone()['alignedSequence'])
	
def db_getUnalSequence(bacterium,gene,allele):
	e = conn.cursor()
	e.execute("SELECT sequence FROM alleles WHERE bacterium = ? AND gene = ? AND alleleVariant = ?",(bacterium,gene,allele))
	return (e.fetchone()['sequence'])
	
def buildConsensus(samFile,chromosomeList,filterScore,max_xM):

	chromosomes = {}
	chromosomesLen = {}

	strr=''
	filterSpec = chromosomeList.keys()[0].split('_')[0].strip()
	
	fili = open(samFile,'r')
	devnull = open('/dev/null', 'w')
	child = subprocess.Popen("samtools view -h - | grep "+filterSpec,shell=True, stdout=subprocess.PIPE, stdin = fili)
	# samFileContentAK = StringIO(child.communicate()[0]) 
	
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
		
			# strr+=line
	
	
	devnull = open('/dev/null', 'w')
	child = subprocess.Popen("samtools view -bS - | samtools sort -o - - | samtools mpileup - ",shell=True, stdout=subprocess.PIPE, stdin = subprocess.PIPE,stderr = devnull)
	out = StringIO(child.communicate(strr)[0]) 
	
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
				#print "EXd",chromosome,nucleotide,chr
				escape = 1
				continue
				
			if escape: 
				if chr.isdigit(): 
					escapeLen = int(chr)+1
					#print "exca",chromosome,nucleotide,chr
				escape = 0 
				
			if escapeLen > 0: 
				escapeLen -=1
				#print "exca Extend",chromosome,nucleotide,chr
				continue
			
			elif chr.upper() in ldict: 
				#print "DO::",chromosome,nucleotide,chr
				ldict[chr.upper()] += 1
			
			
		#if nucleotide <=100 and nucleotide >= 80: continue
		chromosomes[chromosome][nucleotide] = max(ldict, key=ldict.get)
	
	del out
	seqRec=[]

	for chromo,nucleots in chromosomes.items(): 
		sys.stdout.flush()
		
		sequen = ""
		lastSet = 0
		#print ">"+chromo, chromosomesLen[chromo]
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
				rSequen[i] = dbSequen[i]
				cIndex+=1
				
			elif rSequen[i] != dbSequen[i]: 
				rSequen[i] = rSequen[i].lower()
				SNPs+=1
			i+=1
		sequen = ''.join(rSequen)
		
			
		seqRec.append(SeqRecord(Seq(sequen,IUPAC.unambiguous_dna),id=chromo, description = 'CI::'+str(cIndex)+'_SP::'+str(SNPs)))
	print '\r',
	return seqRec
 	
parser = argparse.ArgumentParser()
parser.add_argument("file", help="BAM file containing the sequences")
#parser.add_argument("-i,--text", help="Information on the sample" action="store_true")
parser.add_argument("-o","--organism", help="focus on ORGANISM only (mlstkey)")

parser.add_argument("-p","--penalty", help="penalty for not found allele", default=100, type=int)
parser.add_argument("--minscore", help="minimum score to match (absolute val!)", default=30, type=int)

parser.add_argument("--silent", help="Silent Mode", action='store_true')
parser.add_argument("--out_folder", help="Output Folder")

parser.add_argument("--mincoverage", help="minimum coverage to produce FASTQ assembler files", default=5, type=int)
parser.add_argument("--min_read_len", help="minimum length for bowtie2 reads", default=90, type=int)
parser.add_argument("--max_xM", help="maximum XM", default=20, type=int)
parser.add_argument("--min_accuracy", help="minimum percentage of allele-reconstruction for each organism to pass", default=0.33, type=float)

parser.add_argument("--present_genes", help="percentage of genes needed to be present in organism, default: 100", default=100, type=int)
parser.add_argument("--log", help="generate logs", action="store_true") 
parser.add_argument("-d","--database", help="database file", default="bsb.db")
parser.add_argument("-a","--all_sequences", help="output also known sequences", action="store_true")

args=parser.parse_args()

#PREPARE
conn = sqlite3.connect(args.database)
conn.row_factory = sqlite3.Row
 

fil = open(args.file,'r')
cel = {}
geneMaxiMin = {}
sequenceBank = {}
ignoredReads = 0

c = conn.cursor()
d = conn.cursor()

fileName = (args.file.split('/'))[-1].split('.')[0]

if args.out_folder:
	workUnit = args.out_folder	
else:
	if not os.path.isdir(fileName): os.mkdir(fileName)
	workUnit = fileName	

devnull = open('/dev/null', 'w')
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
		 
		if (args.organism and species == args.organism) or not args.organism:
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


 
print '\r\n\033[94m'+('  '+fileName+'  ').center(70,'-')+'\033[0m'

for speciesKey,species in cel.items():
	
	tVar = dict([(row['geneName'],0) for row in  c.execute("SELECT geneName FROM genes WHERE bacterium = ?",(speciesKey,))])
	
	#GENE PRESENCE 
	if len(tVar) < len(species.keys()):
		print "\033[43m\033[30mDatabase is broken. Exiting\033[0m",len(tVar),len(species.keys())
		sys.exit(0)
	
	for sk in species.keys():
		tVar[sk] = 1
	vals = sum([t for t in tVar.values()])
	#
	
	# NOT ENOUGH GENES
	if int((float(vals)/float(len(tVar)))*100) < args.present_genes:
		print "\033[91m", speciesKey,"\033[0m\t: ", str(vals)+' genes out of '+str(len(tVar))+' MLST targets'
		print "\t\t   Missing: \033[91m"+', '.join([sk for (sk,v) in tVar.items() if v == 0])+'\033[0m'
		sys.stdout.flush()
	# ENOUGH GENES
	else:
		print "\033[92m", speciesKey,"\033[0m\t: ", str(vals)+' genes out of '+str(len(tVar))+' MLST targets'
		print "\t\t   Missing: \033[91m"+', '.join([sk for (sk,v) in tVar.items() if v == 0])+'\033[0m'
		
		print "\r\n  "+"Gene".ljust(7)+"Coverage".rjust(10)+"Score".rjust(6)+"Hits".rjust(5)+" Allele(s)".ljust(40)
		sys.stdout.flush()
		# profileTrack = []
		for geneKey, geneInfo in species.items(): #geni
		
			minValue = max([avg for (val,leng,avg) in geneInfo.values()])
			
			aElements = {}
			tmp=[]
			for k,(val,leng,avg) in geneInfo.items():
				if avg == minValue:
					aElements[k]=(val,leng,avg)
					tmp.append(k)
			tmp = ",".join(tmp)
			
			sequenceKey = speciesKey+'_'+geneKey
			c.execute("SELECT LENGTH(sequence) as L FROM alleles WHERE bacterium = ? AND gene = ? ORDER BY L DESC LIMIT 1", (speciesKey,geneKey))
			genL = c.fetchone()['L']
			#Removed: de-comment to have sequences and quality for FASTQ generation
				# coverage = sum([len(x) for (x,q) in sequenceBank[sequenceKey].values()])
			# print [x for (x,q) in sequenceBank[sequenceKey].values()]
			coverage = sum([x for x in sequenceBank[sequenceKey].values()])
			
			if coverage >= args.mincoverage * genL: 
				color = "\033[92m"
				#Removed: de-comment to have seque nces and quality for FASTQ generation
					# fqfil = open(workUnit+'/'+sequenceKey+'.fasta','a')
					# for sequenceSpec,(sequence,quality) in sequenceBank[sequenceKey].items():
					#	fqfil.write('>'+sequenceSpec+'\r\n')  
					#	fqfil.write(sequence+'\r\n')   
					# fqfil.close()
			else:  color = "\033[93m"
			
			print "  "+color+geneKey.ljust(7)+"\033[0m"+str(round(float(coverage)/float(genL),2)).rjust(10)+"\033[95m"+str(minValue).rjust(6)+str(aElements.itervalues().next()[1]).rjust(5)+"\033[0m"+"\033[94m",tmp.ljust(40)+"\033[0m"
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
		
		print "      Consensous Sequence Build"
			
		l = [[(speciesKey+'_'+g1+'_'+k,db_getUnalSequence(speciesKey,g1,k)) for k,(val,leng,avg) in g2.items() if avg == max([avg1 for (val1,leng1,avg1) in g2.values()])][0] for g1,g2 in species.items()] 
		consenSeq = buildConsensus(args.file, dict(l),args.minscore,args.max_xM)
		
		#_V2
		# l = [[(speciesKey+'_'+g1+'_'+k,db_getUnalSequence(speciesKey,g1,k)) for k,(val,leng,avg) in g2.items() if avg == max([avg1 for (val1,leng1,avg1) in g2.values()])] for g1,g2 in species.items()] 
		# consenSeq = buildConsensus(samFileContent.getvalue().split('\n'), dict(itertools.chain(l)),args.minscore,args.max_xM)
		
		newProfile = 0
		finWrite = 1
		#print "      "+"Gene".ljust(6)+'Ref.'.ljust(7)+"Length".rjust(10)+"Ns".rjust(10)+"SNPs".rjust(10)
		print "\r\n  "+"Gene".ljust(6)+"Ref.".ljust(7)+"Length".rjust(7)+"Ns".rjust(7)+"SNPs".rjust(10)+"Accuracy".rjust(9)+"Notes".rjust(10)
		for l in consenSeq:
			holes = str(l.description.split('_')[0].split('::')[1])
			snps = int(l.description.split('_')[1].split('::')[1])
			leng = str(len(l.seq))
			leng_ns = str(round(1-float(holes)/float(leng),4)*100)+' %'
			l.seqLen = len(l.seq)
			
			# If there's a gene with low accuracy, the whole organism is discarded for this sample
			if (1-float(holes)/float(leng)) <= args.min_accuracy: finWrite = 0
			
			if snps > 0:
				seqFind = sequenceFind(speciesKey,l.seq)
				if seqFind: newAllele = seqFind
				else:
					newAllele = 'NEW'
					newProfile = 1
					
			else:
				newAllele = '--'
				if not args.all_sequences: l.seq = ''
			
			print "  "+(l.id.split('_')[1]).ljust(6)+(l.id.split('_')[2]).ljust(7)+leng.rjust(7)+holes.rjust(7)+str(snps).rjust(10)+leng_ns.rjust(9)+newAllele.rjust(10)
			
		
		print ''
		 
		
		# Output File Track:
		# bacterium_name	KP	MLST_PROFILE_CODE	gene::sequence	gene::sequence...
		
		# for record in consenSeq:
			# if int(record.description.split('_')[1].split('::')[1]) > 0: 	#DIFFERENT SEQUENCE (SNPs > 0)
				# if not sequenceExists(speciesKey,record.seq): 						#NEW PROFILE MARKER 
					# newProfile=1
					
		if finWrite:
			print '  Reconstruction Successful > Write'
			profil = open(workUnit+'/'+fileName+'.txt','a')	 
			#recd.description.split('_')[0] = "CI::Confidence_Index_value << Ns"
			#recd.description.split('_')[1] = "SN::SNPs_values << SNPs"
			profil.write(speciesKey+'\t'+fileName+'\t'+"\t".join([recd.id+"::"+str(recd.seq)+'::'+str(round(1-float(recd.description.split('_')[0].split('::')[1])/float(recd.seqLen),4)*100)+'::'+str(round(float(recd.description.split('_')[1].split('::')[1])/float(recd.seqLen),4)*100) for recd in consenSeq])+'\r\n')
			
			profil.close()
		else: print '  Reconstruction Failed  ('+str(round(args.min_accuracy*100,2))+'% min.accuracy) > \033[91mSkip\033[0m'
		print ''
		
	del cel[speciesKey]
	gc.collect()
			
		# else: 
			# profLook = defineProfile([rec.id for rec in consenSeq])[0]
			# if profLook[1] == 100: # KP
				# profil.write(speciesKey+"\tKP\t"+str(profLook[0])+"\t"+"\t".join([recd.id+"::"+str(recd.seq) for recd in consenSeq])+'\r\n')
			# else:  
				# profil.write(speciesKey+"\tNP\t"+str(profLook[0])+'::'+str(profLook[1])+"\t"+"\t".join([recd.id+"::"+str(recd.seq) for recd in consenSeq])+'\r\n')
				
		
		
		
		#print defineProfile(["haemophilus_mdh30","haemophilus_frdB152","haemophilus_recA31","haemophilus_pgi1","haemophilus_atpG50","haemophilus_adk170"])

		## GLOBAL SEQUENCE 
		# for geneKey, geneInfo in sorted(species.items(), key= lambda x: x[0]): #for each (ordered) gene of the organism
	
			# alleles = [k for k,(val,leng,avg) in geneInfo.items() if avg == min([avg1 for (val1,leng1,avg1) in geneInfo.values()])]
			# kt=0
			
			# for allele in alleles:
				# kt = kt+1
				# if alleleInManyProfile(speciesKey,geneKey,allele,matchingProfiles): 
					# print ' \033[92m',(str(geneKey)+str(allele)).ljust(10),'\033[0m\tFOUND! Adding...'
					# sampleSequence = sampleSequence + db_getSequence(speciesKey,geneKey,allele)
					# break
				# else:
					# if kt == len(alleles):
						# sampleSequence = sampleSequence + db_getSequence(speciesKey,geneKey,allele)
						# print ' \033[94m',(str(geneKey)+str(allele)).ljust(10),'\033[0m\tFOUND! Adding...'
					# else: 
						# print ' \033[91m',(str(geneKey)+str(allele)).ljust(10),'\033[0m\tNOT FOUND'
		
		# outSequence = []
		# outSequence.append(SeqRecord(Seq(sampleSequence, IUPAC.unambiguous_dna), id = 'sample_'+str(speciesKey)+' [closest profile(s): '+repr(matchingProfiles)+']', description = speciesKey))
		# SeqIO.write(outSequence, 'seq_'+speciesKey+'.fasta', "fasta")
			 
conn.close() 

if len(cel): print '\033[92m'+'[ - Completed - ]'.rjust(70,' ')+'\033[0m'
#print cel
