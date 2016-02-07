import sys,os,subprocess,sqlite3,argparse,difflib,math,itertools
from StringIO import StringIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Align.Applications import MuscleCommandline

 
class bcolors:
	HEADER = '\033[1;95m'
	OKBLUE = '\033[1;94m'
	OKGREEN = '\033[1;92m'
	WARNING = '\033[1;93m'
	RED = '\033[1;91m'
	CYAN = '\033[0;37m'
	ENDC = '\033[0m'
	OKGREEN2 = '\033[42m\033[30m' 


	

parser = argparse.ArgumentParser()
parser.add_argument("folder", help="folder containing .txt sametagenome files")
parser.add_argument("-d","--database", help="database", default="bsb.db")
parser.add_argument("--meta", help="metadata file (CSV)")
parser.add_argument("--idField", help="field containig the 'sampleID' value", default=0, type=int)
parser.add_argument("-z","--allele_max_snps", help="edit distance threshold to pass to call a new MLST allele", default=100, type=int)
 
parser.add_argument("--outseqformat", choices=['A', 'A+', 'B', 'B+', 'C','C+'], help="A  : Concatenated Fasta (Detected STs)\n \
		    A+ : Concatenated Fasta (All STs)\n \
		    B  : Single loci (New loci only)\n \
		    B+ : Single loci (All loci)\n \
		    C  : CSV ST Table")
parser.add_argument('--filter', help="list of subset of organisms to focus on (comma separated)")
parser.add_argument("-j", help="output sequences with a list of metadata (comma separated)")
parser.add_argument("--jgroup", help="output sequences, group by ST", action="store_true")

parser.add_argument("--muscle_path", default="tools/muscle")

args=parser.parse_args()

try:
	conn = sqlite3.connect(args.database)
except IOError: 
	print "IOError: unable to access "+args.database+"!"
conn.row_factory = sqlite3.Row
cursor = conn.cursor()

def defineProfile(geneList):
	
	recs=[]
	for allele in geneList:
		e = conn.cursor()
		e.execute("SELECT recID FROM alleles WHERE bacterium||'_'||gene||'_'||alleleVariant = ?",(allele,))
		result = e.fetchone()
		
		if result:
			recs.append(str(result['recID']))

	return [(row['profileCode'],int((float(row['T'])/float(len(recs)))*100)) for row in e.execute("SELECT profileCode, COUNT(*) as T FROM profiles WHERE alleleCode IN ("+','.join(recs)+") GROUP BY profileCode HAVING T = (SELECT COUNT(*)  FROM profiles WHERE alleleCode IN ("+','.join(recs)+") GROUP BY profileCode ORDER BY COUNT(*) DESC LIMIT 1) ORDER BY T DESC")] if result else [(0,0)]
	
def sequenceExists(bacterium,sequence):
	
	e = conn.cursor()
	e.execute("SELECT 1 FROM alleles WHERE sequence = ? AND bacterium = ?",(str(sequence),bacterium))
	return (len(e.fetchall()) > 0)

def sequenceLocate(bacterium,sequence):
	
	e = conn.cursor()
	e.execute("SELECT alleleVariant FROM alleles WHERE sequence = ? AND bacterium = ?",(str(sequence),bacterium))
	return str(e.fetchone()['alleleVariant'])

def sequencesGetAll(bacterium,gene):
	
	e = conn.cursor()
	e.execute("SELECT sequence,alleleVariant FROM alleles WHERE gene = ? AND bacterium = ?",(gene,bacterium))
	return dict((x['alleleVariant'],x['sequence']) for x in e.fetchall())

def stringDiff(s1,s2):
	c =0 
	for a,b in zip(s1,s2):
		if a!=b: c+=1
	return c
		
cel = {}

if not os.path.isdir(args.folder+'/merged'): os.makedirs(args.folder+'/merged')

 
for file in os.listdir(args.folder):
	
	if file.split('.')[-1] != 'nfo': continue 
 
	
	for line in open(args.folder+'/'+file,'r'):
		organism = line.split()[0] 
		sampleName = line.split()[1] 
		genes =line.split()[2::] 
		
		#apply the species filter
		if args.filter and organism not in args.filter.split(','): continue
		
		if organism not in cel: cel[organism] = []
		cel[organism].append((dict((x.split('::')[0],(x.split('::')[1].upper(),x.split('::')[2],x.split('::')[3])) for x in genes),sampleName))

for bacterium,bactRecord in cel.items(): #For each bacterium:

	#bactRecord = ( {gene_allele: (sequence,accuracy,snps) , ...} , sampleName ) 
	
	print '\r\n----------------------------------------------------------------------\r\n'+bacterium+'\r\n----------------------------------------------------------------------\r\n'
	 
	profil = open(args.folder+'/merged/'+bacterium+'_ST.txt','w')
	
	#profil.write(bacterium+'\r\n')
	
	stringBase = {}
	
	oldProfiles = {}
	genesBase={}
	profilesBase={}
	encounteredProfiles={}
	isolates=[]
	newSequences = {} #contains the new seqs (reconstructed) for further alignment
	mainGeneStructure = {} #key:label,value:sequence. 
	
	
	cursor.execute("SELECT profileCode FROM profiles WHERE bacterium = ? ORDER BY profileCode DESC LIMIT 1",(bacterium,))

	# lastProfile = int(cursor.fetchone()['profileCode'])
	lastProfile = 100000
	# lastGenes = dict((row['gene'],row['maxGene']) for row in cursor.execute("SELECT gene, MAX(alleleVariant) as maxGene FROM alleles WHERE bacterium = ? GROUP BY gene",(bacterium,)))
	lastGenes = dict((row['gene'],100000) for row in cursor.execute("SELECT gene, MAX(alleleVariant) as maxGene FROM alleles WHERE bacterium = ? GROUP BY gene",(bacterium,)))
	
	#Get all KNOWN profiles for that bacterium
	
	for row in cursor.execute("SELECT profileCode,gene,alleleVariant FROM profiles,alleles WHERE alleleCode = alleles.recID AND alleles.bacterium = ?",(bacterium,)):
		if row['profileCode'] not in oldProfiles: oldProfiles[row['profileCode']] = [0,{}]
		oldProfiles[row['profileCode']][1][row['gene']] = row['alleleVariant']
	
	for bacteriumLine,sampleRecord in bactRecord: #for each entry of that bacterium:
		 
		#bacteriumLine is a dict. Key: genes, Values: sequences. No sequence = no news.
		#>> {gene_allele: (sequence,accuracy) , ...} 
		profileLine = {}
		newAlleles = []
		flagRecurrent = False
		sum_of_accuracies = 0.0
		
		for geneLabel,(geneSeq,geneAccur,percent_snps) in bacteriumLine.items(): #for each gene of the entry
			 #TODO: change seq_recog
			geneOrganism,geneName,geneAllele = geneLabel.split('_')
			sum_of_accuracies += float(geneAccur)
			if geneSeq == '' or sequenceExists(bacterium,geneSeq):
				# WE HAVE A DATABASE SEQUENCE
				#print "WE HAVE A DATABASE SEQUENCE"
				#geneSeq = cursor.execute("SELECT sequence")
				if geneSeq != '': geneAllele = sequenceLocate(bacterium,geneSeq)
				profileLine[geneName] = (geneAllele,0)
				
			elif geneSeq in genesBase:
				profileLine[geneName] = (genesBase[geneSeq].split('_')[2],2)
				flagRecurrent = True #a new allele is recurring
 
			elif geneSeq not in genesBase:
				## WE HAVE A NEW SEQUENCE
				 
				
				geneCategoryCode = 1 #default: blue (new allele, accepted)
				
				if args.allele_max_snps != None: 
					geneCategoryCode = 3 #default becomes now not accepted
					
					for refCode,refSeq in sequencesGetAll(bacterium,geneName).items():
						if stringDiff(geneSeq,refSeq) <= args.allele_max_snps:
						    #print geneName,refCode,stringDiff(geneSeq,refSeq)
						    geneCategoryCode = 1 # if match allele_max_snps: accept
						    break
						    
					
				geneNewAlleleNumber = str(lastGenes[geneName]+1)
				lastGenes[geneName]+=1
				
				geneNewLabel = geneOrganism+'_'+geneName+'_'+geneNewAlleleNumber
				genesBase[geneSeq] = geneNewLabel 
					#print "\t New Sequence for "+geneName+" -> "+geneNewLabel
				profileLine[geneName] = (geneNewAlleleNumber,geneCategoryCode) # 1: blue (new allele, accepted)
											       # 3: red  (new allele, not accepted)
				newAlleles.append(geneName)
				
					
				if geneName not in newSequences: newSequences[geneName] = []
				newSequences[geneName].append(SeqRecord(Seq(geneSeq,IUPAC.unambiguous_dna),id=geneNewLabel, description = ''))
 
		
		meanAccuracy = sum_of_accuracies / float(len(bacteriumLine))
		if len(newAlleles) == 0: 
		## Existent MLST profile -> Look for it (--> ISOLATES)
		
			#Tries to define an existing MLST profile with the alleles
			if not flagRecurrent:
				tryDefine = defineProfile([bacterium+'_'+k+'_'+v[0] for k,v in profileLine.items()])
				#lopo
				if tryDefine and tryDefine[0][1] == 100: 
			
					#encounteredProfiles[tryDefine[0][0]] = [profileLine,1,0]
					oldProfiles[tryDefine[0][0]][0]+=1
					isolates.append((tryDefine[0][0],meanAccuracy,sampleRecord)) 
					continue
			
			
			foundExistant = 0
			for key,(element,abundance,isNewProfile) in encounteredProfiles.items():
				if [k+str(v[0]) for k,v in sorted(profileLine.items())] == [k+str(v[0]) for k,v in sorted(element.items())]:
					foundExistant = key
					
			if foundExistant:
				encounteredProfiles[foundExistant][1] += 1
				isolates.append((foundExistant,meanAccuracy,sampleRecord)) 
			else:
				lastProfile+=1
				encounteredProfiles[lastProfile] = [profileLine,1,2]
				isolates.append((lastProfile,meanAccuracy,sampleRecord)) 
		else:
			# THIS IS A NEW PROFILE
			lastProfile+=1
			
			profileCategoryCode = 1 
			if args.allele_max_snps != None: 
				for k,(v,cat) in profileLine.items():
				    if cat == 3: #if rejectable allele
					profileCategoryCode = 3 #rejectable profile
					break
				      
				      
			
			encounteredProfiles[lastProfile] = [profileLine,1,profileCategoryCode]
			if profileCategoryCode != 3: isolates.append((lastProfile,meanAccuracy,sampleRecord)) 
		
		
			
		
		# PROFILE LINE: dictionary['gene']: allele, hits, color
		#	-- ExistingCode | Effect --
		#           1       |  GREEN  -> New, with new alleles
		#           2       |  YELLOW -> New, with old allelese
		#           3	    |  RED    -> New, with new alleles some of which rejectable	
		#                   | 
	#Old profiles
	
	
	
	
	profil.write('ST\t'+'\t'.join([x for x in sorted(lastGenes.keys())] )+'\r\n')
	print 'KNOWN MLST profiles found:\nST\t'+'\t'.join([x for x in sorted(lastGenes.keys())])+'\tHits'
	
	# OLD PROFILES
	
	for profileCode,(hits,profile) in oldProfiles.items():
		profil.write(str(profileCode)+'\t'+'\t'.join([str(v) for k,v in sorted(profile.items())])+'\r\n')
		if hits > 0: 
			print bcolors.RED+str(profileCode)+bcolors.ENDC + '\t' + '\t'.join([str(v) for k,v in sorted(profile.items())])+'\t'+str(hits)
			sys.stdout.flush()
			
	#NEW PROFILES (ACCEPTED)
	
	print '\n\nNEW MLST profiles found:\nST\t'+'\t'.join([x for x in sorted(lastGenes.keys())])+'\tHits'	
	for profileID,(profile,hits,profileCategoryCode) in encounteredProfiles.items():
		
		if profileCategoryCode not in [1,2]: continue
	      
		if profileCategoryCode == 1: profileNumber = bcolors.OKGREEN+str(profileID)+bcolors.ENDC
		elif profileCategoryCode == 2: profileNumber = bcolors.WARNING+str(profileID)+bcolors.ENDC
		
		print profileNumber + '\t' + '\t'.join([bcolors.OKBLUE+str(v[0])+bcolors.ENDC if v[1] == 1 else bcolors.OKGREEN+str(v[0])+bcolors.ENDC if v[1] == 2 else bcolors.RED+str(v[0])+bcolors.ENDC if v[1] == 3 else str(v[0]) for k,v in sorted(profile.items())]) + '\t' + str(hits)
		sys.stdout.flush()
			
		profil.write(str(profileID)+'\t'+'\t'.join([str(v[0]) for k,v in sorted(profile.items())])+'\n')

	#NEW PROFILES (REJECTED)
	
	print '\n\nREJECTED NEW MLST profiles, SNPs threshold: '+str(args.allele_max_snps)+'\nST\t'+'\t'.join([x for x in sorted(lastGenes.keys())])+'\tHits'
	for profileID,(profile,hits,profileCategoryCode) in encounteredProfiles.items():
		
		if profileCategoryCode not in [3]: continue
	        profileNumber = bcolors.CYAN+str(profileID)+bcolors.ENDC
	        
		print str(profileID) + '\t' + '\t'.join([bcolors.OKBLUE+str(v[0])+bcolors.ENDC if v[1] == 1 else bcolors.OKGREEN+str(v[0])+bcolors.ENDC if v[1] == 2 else bcolors.RED+str(v[0])+bcolors.ENDC if v[1] == 3 else str(v[0]) for k,v in sorted(profile.items())]) + '\t' + str(hits)
		sys.stdout.flush()
			
		#profil.write(str(profileID)+'\t'+'\t'.join([str(v[0]) for k,v in sorted(profile.items())])+'\n')
	
	profil.close()
	
	
	#ISOLATES FILE OUTPUT 
	isolafil = open(args.folder+'/merged/'+bacterium+'_isolates.txt','w') 
	identifiers = {}
	p1line = 0
	keys=[]
	
	if args.meta:
		for line in open(args.meta):
		      if line == '': continue
		      if not p1line: 
			      p1line=1
			      keys = [str(x).strip() for x in line.split('\t')]
		      else:
			      l = line.strip().split('\t') 
			      if len(l) == len(keys):
				  identifiers[l[args.idField]] = dict((keys[i],l[i]) for i in range(0,len(keys)))  
				  #print "put",l[args.idField]
	  
	
	isolafil.write('ST\tBreadth of coverage\t'+'\t'.join(keys)+'\n')
	
	
	
	STmapper={} #used to keep track of STs to output in form of concatenated sequences
	
	for profileST,meanAccur,sampleName in isolates:
		
		if profileST not in STmapper: STmapper[profileST] = [] 
		if sampleName.endswith('.fna'): sampleName = sampleName.split('.')[0]
		
		if sampleName in identifiers:  
			strl=[]
			for ky in keys:
				strl.append(identifiers[sampleName][ky])
			isolafil.write(str(profileST)+'\t'+str(round(meanAccur,2))+'\t'+'\t'.join(strl)+'\n')
			
			STmapper[profileST].append(identifiers[sampleName])
		else: 
			isolafil.write(str(profileST)+'\t'+str(round(meanAccur,2))+'\t'+str(sampleName)+'\n')
			
			#if args.j: STmapper[profileST] = dict((virtKey,'-') for virtKey in args.j)
			STmapper[profileST].append({'sampleID':sampleName})
			
			
	isolafil.close()
		
	
	#SEQUENCES OUTPUT
	if args.outseqformat:
	
		if args.outseqformat == 'B':
			SeqIO.write(sorted( list(itertools.chain(*newSequences.values()))  ,key=lambda x: x.id),args.folder+'/merged/'+bacterium+'_sequences.fna', "fasta")
		
		seqTable={}
		preaLignTable={}
		for row in cursor.execute("SELECT gene,alleleVariant,sequence FROM alleles WHERE bacterium = ? ORDER BY bacterium,gene,alleleVariant",(bacterium,)):
			
			label = bacterium+'_'+row['gene']+'_'+str(row['alleleVariant'])
			if row['gene'] not in preaLignTable: preaLignTable[row['gene']] = []
			preaLignTable[row['gene']].append(SeqRecord(Seq(row['sequence'], IUPAC.unambiguous_dna), id = label, description=''))

		for seqGene,seqList in newSequences.items():
			if seqGene not in preaLignTable: preaLignTable[seqGene] = []
			for seqElement in seqList:
				preaLignTable[seqGene].append(seqElement)
					
		if args.outseqformat == 'B+':
			SeqIO.write(sorted( list(itertools.chain(*preaLignTable.values()))  ,key=lambda x: x.id),args.folder+'/merged/'+bacterium+'_sequences.fna', "fasta")	
			
				
		if args.outseqformat == 'C':
			
			seqfile = open(args.folder+'/merged/'+bacterium+'_sequences.txt','w')
			
			
			nalign_Table= dict((k.id,k.seq) for k in list(itertools.chain(*preaLignTable.values()))) 
			
			  
			seqfile.write('ST\t'+'\t'.join([str(x) for x in sorted(lastGenes.keys())] )+'\r\n')
			for profileCode,(hits,profile) in oldProfiles.items(): #for each old profile
				if hits>0 or args.outseqformat == 'C+':
					seqfile.write(str(profileCode)+'\t'+'\t'.join([str(nalign_Table[bacterium+'_'+gen+'_'+str(alle)]) for gen,alle in sorted(profile.items())] )+'\r\n')

			for profileCode,(profile,hits,isNewProfile) in encounteredProfiles.items(): #for each encountered profile
				if isNewProfile == 3: continue #rejected profiles
			      
				seqfile.write(str(profileCode)+'\t'+'\t'.join([str(nalign_Table[bacterium+'_'+gen+'_'+str(alle[0])]) for gen,alle in sorted(profile.items())] )+'\r\n')

			seqfile.close()
 
		if args.outseqformat in ['A','A+']: #sequences, merged

			for gene,seqs in preaLignTable.items(): #for each gene
				
				
				cS = StringIO()
				tld=[]
				for seq in seqs:
					if len(seq) not in tld: tld.append(len(seq))
				
				if len(tld) > 1: #more than one length: need to align!
				
					print(('\tAligning Sequences ['+gene+']:').ljust(50)),
					sys.stdout.flush()
					SeqIO.write(seqs,cS, "fasta")
					muscle_cline = MuscleCommandline(args.muscle_path)
					stdout, stderr = muscle_cline(stdin=cS.getvalue())
					for sequence in SeqIO.parse(StringIO(stdout), "fasta"):
						seqTable[sequence.id] = str(sequence.seq)
				else:
					print(('\tAligned Sequences ['+gene+'] ('+str(len(tld))+'):').ljust(50)),
					for seq in seqs:
						seqTable[seq.id] = str(seq.seq) 
						
				print(bcolors.OKGREEN+'[ - Done - ]'+bcolors.ENDC)
				sys.stdout.flush()
	
			phyloSeq = []
			#OLD 
			for profileCode,(hits,profile) in oldProfiles.items(): #for each old profile
				stSeq = ''
				 

				if hits>0: #old, detected, present in the samples

					for gen,all in sorted(profile.items()):
						#print gen,all 
						stSeq+=str(seqTable[bacterium+'_'+gen+'_'+str(all)])
					
					if args.j:  
						listofkeys = dict((k,[]) for k in args.j.split(','))
						
						if profileCode in STmapper: 
							prog = 0
							for i in [x for x in STmapper[profileCode]]:
								if args.jgroup:
									descriptionString='n='+str(hits)
									for (kl,v) in i.items():
										if kl in listofkeys.keys(): listofkeys[kl].append(v)
									descriptionString+= ''.join([kll+'{'+'|'.join(ell)+'}' for kll,ell in listofkeys.items()])
								else:
									prog+=1
									descriptionString = '-'.join([kll+'{'+str(ell)+'}' for kll,ell in i.items()  if kll in args.j.split(',')])
									phyloSeq.append(SeqRecord(Seq(stSeq,IUPAC.unambiguous_dna),id=bacterium+'_ST'+str(profileCode)+'_'+str(prog)+'_'+descriptionString, description = descriptionString))
									
						#output one sequence, for the current ST with all the metadata grouped
						if args.jgroup: phyloSeq.append(SeqRecord(Seq(stSeq,IUPAC.unambiguous_dna),id=bacterium+'_ST'+str(profileCode)+'_'+descriptionString, description = descriptionString))
					else:
						for profileInstance in STmapper[profileCode]:
							phyloSeq.append(SeqRecord(Seq(stSeq,IUPAC.unambiguous_dna),id=profileInstance['sampleID'], description = ''))
				
				elif args.outseqformat == 'A+': #old, non present, but required to be added
					for gen,all in sorted(profile.items()):
						#print gen,all 
						stSeq+=str(seqTable[bacterium+'_'+gen+'_'+str(all)])
					phyloSeq.append(SeqRecord(Seq(stSeq,IUPAC.unambiguous_dna),id="ST_"+str(profileCode), description = ''))

						
			#NEW
			for profileCode,(profile,hits,isNewProfile) in encounteredProfiles.items(): #for each encountered profile
				
				if isNewProfile == 3: continue #rejected profiles
				
				stSeq = ''
				for gen,all in sorted(profile.items()):
					stSeq+=str(seqTable[bacterium+'_'+gen+'_'+all[0]])
					
				#presence of metadata
				if args.j:
					listofkeys = dict((k,[]) for k in args.j.split(','))
						
					if profileCode in STmapper: 
						prog = 0
						for i in [x for x in STmapper[profileCode]]:
							if args.jgroup:
								descriptionString='n='+str(hits)
								for (kl,v) in i.items():
									if kl in listofkeys.keys(): listofkeys[kl].append(v)
								descriptionString+= ''.join([kll+'{'+'|'.join(ell)+'}' for kll,ell in listofkeys.items()])
							else:
								#for (kl,v) in i.items():
								prog+=1
								descriptionString = '-'.join([kll+'{'+str(ell)+'}' for kll,ell in i.items() if kll in args.j.split(',')])
								phyloSeq.append(SeqRecord(Seq(stSeq,IUPAC.unambiguous_dna),id=bacterium+'_ST'+str(profileCode)+'_'+str(prog)+'_'+descriptionString, description = descriptionString))
					#output one sequence, for the current ST with all the metadata grouped
					if args.jgroup: phyloSeq.append(SeqRecord(Seq(stSeq,IUPAC.unambiguous_dna),id=bacterium+'_ST'+str(profileCode)+'_'+descriptionString, description = descriptionString))
				else:
					for profileInstance in STmapper[profileCode]: 
						phyloSeq.append(SeqRecord(Seq(stSeq,IUPAC.unambiguous_dna),id=profileInstance['sampleID'], description = ''))
			
			SeqIO.write(phyloSeq,args.folder+'/merged/'+bacterium+'_sequences.fna', "fasta")

print "Completed"
