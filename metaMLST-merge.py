#!/usr/bin/env python
try: 
	import math
	import sqlite3
	import subprocess
	import time
	import itertools
	import gc
	import os
	import sys
	from metaMLST_functions import *
	from StringIO import StringIO
except ImportError as e:
	print "Error while importing python modules! Remember that this script requires: sys,os,subprocess,sqlite3,argparse,re"
	sys.exit(1)
 

try:
	from Bio import SeqIO
	from Bio.Seq import Seq
	from Bio.SeqRecord import SeqRecord
	from Bio.Alphabet import IUPAC
	from Bio.Align.Applications import MuscleCommandline
except ImportError as e:
	metamlst_print("Failed in importing Biopython. Please check Biopython is installed properly on your system!",'FAIL',bcolors.FAIL)
	sys.exit(1)


 
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
		description='Detects the MLST profiles from a collection of intermediate files from MetaMLST.py')

parser.add_argument("folder", help="Path to the folder containing .nfo MetaMLST.py files")
parser.add_argument("-d", metavar="DB_PATH", help="MetaMLST SQLite Database File (created with metaMLST-index)", required=True)
parser.add_argument("--filter", metavar="species1,species2...", help="Filter for specific set of organisms only (METAMLST-KEYs, comma separated. Use metaMLST-index.py --listspecies to get MLST keys)")
parser.add_argument("-z", metavar="ED", help="Maximum Edit Distance from the closest reference to call a new MLST allele. Default: 5", default=5, type=int)

parser.add_argument("--meta", metavar="METADATA_PATH", help="Metadata file (CSV)")
parser.add_argument("--idField", help="Field number pointing to the 'sampleID' value in the metadata file", default=0, type=int) 
parser.add_argument("--outseqformat", choices=['A', 'A+', 'B', 'B+', 'C','C+'], help="A  : Concatenated Fasta (Only Detected STs)\r\n\
A+ : Concatenated Fasta (All STs)\r\n\
B  : Single loci (Only New Loci)\r\n\
B+ : Single loci (All loci)\r\n\
C  : CSV STs Table [default]")
parser.add_argument("-j", metavar="subjectID,diet,age...", help="Embed a LIST of metadata in the the output sequences (A or A+ outseqformat modes). Requires a comma separated list of field names from the metadata file specified with --meta")
parser.add_argument("--jgroup", help="Group the output sequences (A or A+ outseqformat modes) by ST, rather than by sample. Requires -j", action="store_true")

args=parser.parse_args()

try:
	conn = sqlite3.connect(args.d)
except IOError: 
	print "IOError: unable to access "+args.d+"!"
conn.row_factory = sqlite3.Row
cursor = conn.cursor()


		
cel = {}

if not os.path.isdir(args.folder+'/merged'): os.makedirs(args.folder+'/merged')

#print defineProfile(conn,['ecoli_adk_10','ecoli_fumC_11','ecoli_gyrB_4','ecoli_icd_8','ecoli_mdh_8','ecoli_purA_8','ecoli_recA_2'])
#print defineProfile(conn,['ecoli_adk_21','ecoli_fumC_35','ecoli_gyrB_27','ecoli_icd_6','ecoli_mdh_5','ecoli_purA_5','ecoli_recA_4'])

#sys.exit(0)
 
for file in os.listdir(args.folder):
	
	if file.split('.')[-1] != 'nfo': continue 
 
	
	for line in open(args.folder+'/'+file,'r'):
		organism = line.split()[0] 
		sampleName = line.split()[1] 
		genes =line.split()[2::] 
		
		#apply the species filter
		if args.filter and organism not in args.filter: continue
		
		if organism not in cel: cel[organism] = []
		cel[organism].append((dict((x.split('::')[0],(x.split('::')[1].upper(),x.split('::')[2],x.split('::')[3])) for x in genes),sampleName))

for bacterium,bactRecord in cel.items(): #For each bacterium:


	print bcolors.OKBLUE+('-'*80)+bcolors.ENDC 
	print bcolors.OKBLUE+'|'+bcolors.ENDC+str(db_getOrganisms(conn,bacterium)).center(78)+bcolors.OKBLUE+'|'+bcolors.ENDC
	print bcolors.OKBLUE+('-'*80)+bcolors.ENDC 
	 
	profil = open(args.folder+'/merged/'+bacterium+'_ST.txt','w')
	
	stringBase = {}
	oldProfiles = {}
	genesBase={}
	profilesBase={}
	encounteredProfiles={}
	isolates=[]
	newSequences = {} #contains the new seqs (reconstructed) for further alignment
	mainGeneStructure = {} #key:label,value:sequence. 
	metadataJoinField = 'sampleID'

	cursor.execute("SELECT profileCode FROM profiles WHERE bacterium = ? ORDER BY profileCode DESC LIMIT 1",(bacterium,))

	
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
			if geneSeq == '' or sequenceExists(conn,bacterium,geneSeq):
				# WE HAVE A DATABASE SEQUENCE
				#print "WE HAVE A DATABASE SEQUENCE"
				#geneSeq = cursor.execute("SELECT sequence")
				if geneSeq != '': geneAllele = sequenceLocate(conn,bacterium,geneSeq)
				profileLine[geneName] = (geneAllele,0)
				
			elif geneSeq in genesBase:
				profileLine[geneName] = (genesBase[geneSeq].split('_')[2],2)
				flagRecurrent = True #a new allele is recurring
 
			elif geneSeq not in genesBase:
				## WE HAVE A NEW SEQUENCE
				 
				
				geneCategoryCode = 1 #default: blue (new allele, accepted)
				
				if args.z != None: 
					geneCategoryCode = 3 #default becomes now not accepted
					
					for refCode,refSeq in sequencesGetAll(conn,bacterium,geneName).items():
						if stringDiff(geneSeq,refSeq) <= args.z:
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
				tryDefine = defineProfile(conn,[bacterium+'_'+k+'_'+v[0] for k,v in profileLine.items()])
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
			if args.z != None: 
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
			print bcolors.FAIL+str(profileCode)+bcolors.ENDC + '\t' + '\t'.join([str(v) for k,v in sorted(profile.items())])+'\t'+str(hits)
			sys.stdout.flush()
			
	#NEW PROFILES (ACCEPTED)
	
	print '\n\nNEW MLST profiles found:\nST\t'+'\t'.join([x for x in sorted(lastGenes.keys())])+'\tHits'	
	for profileID,(profile,hits,profileCategoryCode) in encounteredProfiles.items():
		
		if profileCategoryCode not in [1,2]: continue
	      
		if profileCategoryCode == 1: profileNumber = bcolors.WARNING+str(profileID)+bcolors.ENDC
		elif profileCategoryCode == 2: profileNumber = bcolors.OKGREEN+str(profileID)+bcolors.ENDC
		
		print profileNumber + '\t' + '\t'.join([bcolors.OKBLUE+str(v[0])+bcolors.ENDC if v[1] == 1 else bcolors.OKGREEN+str(v[0])+bcolors.ENDC if v[1] == 2 else bcolors.FAIL+str(v[0])+bcolors.ENDC if v[1] == 3 else str(v[0]) for k,v in sorted(profile.items())]) + '\t' + str(hits)
		sys.stdout.flush()
			
		profil.write(str(profileID)+'\t'+'\t'.join([str(v[0]) for k,v in sorted(profile.items())])+'\n')

	#NEW PROFILES (REJECTED)
	
	print '\n\nREJECTED NEW MLST profiles, SNPs threshold (-z): '+str(args.z)+'\nST\t'+'\t'.join([x for x in sorted(lastGenes.keys())])+'\tHits'
	for profileID,(profile,hits,profileCategoryCode) in encounteredProfiles.items():
		
		if profileCategoryCode not in [3]: continue
	        profileNumber = bcolors.OKBLUE+str(profileID)+bcolors.ENDC
	        
		print str(profileID) + '\t' + '\t'.join([bcolors.OKBLUE+str(v[0])+bcolors.ENDC if v[1] == 1 else bcolors.OKGREEN+str(v[0])+bcolors.ENDC if v[1] == 2 else bcolors.FAIL+str(v[0])+bcolors.ENDC if v[1] == 3 else str(v[0]) for k,v in sorted(profile.items())]) + '\t' + str(hits)
		sys.stdout.flush()
			
		#profil.write(str(profileID)+'\t'+'\t'.join([str(v[0]) for k,v in sorted(profile.items())])+'\n')
	
	profil.close()
	print ""
 
	metamlst_print("Outputing results",'...',bcolors.ENDC)
	
	#ISOLATES FILE OUTPUT 
	isolafil = open(args.folder+'/merged/'+bacterium+'_report.txt','w') 
	identifiers = {}
	p1line = 0
	keys=[]

	
	if args.meta:

		for line in open(args.meta):
		      if line == '': continue
		      if not p1line: 
			      p1line=1
			      keys = [str(x).strip() for x in line.split('\t')]
			      metadataJoinField = keys[args.idField]
		      else:
			      l = line.strip().split('\t') 
			      if len(l) == len(keys):

				  identifiers[l[args.idField]] = dict((keys[i],l[i]) for i in range(0,len(keys)))  
				  
	
	isolafil.write('ST\tConfidence\t'+'\t'.join(keys)+'\n')
	
	
	
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
				
					#print(('\tAligning Sequences ['+gene+']:').ljust(50)),
					metamlst_print('Sequences of ['+gene+'] need to be aligned','...',bcolors.ENDC)
					sys.stdout.flush()
					SeqIO.write(seqs,cS, "fasta")
					muscle_cline = MuscleCommandline('muscle')
					stdout, stderr = muscle_cline(stdin=cS.getvalue())
					for sequence in SeqIO.parse(StringIO(stdout), "fasta"):
						seqTable[sequence.id] = str(sequence.seq)
					metamlst_print('Sequences of ['+gene+'] need to be aligned','DONE',bcolors.OKGREEN)
				else:
					#print(('\tAligned Sequences ['+gene+'] ('+str(len(tld))+'):').ljust(50)),
					metamlst_print('Sequences of ['+gene+'] are aligned','OK',bcolors.OKGREEN)
					for seq in seqs:
						seqTable[seq.id] = str(seq.seq) 
						
				#print(bcolors.OKGREEN+'[ - Done - ]'+bcolors.ENDC)
			metamlst_print('Sequences Alignment Completed','DONE',bcolors.OKGREEN)
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
									phyloSeq.append(SeqRecord(Seq(stSeq,IUPAC.unambiguous_dna),id=bacterium+'_ST'+str(profileCode)+'_'+str(prog)+'_'+descriptionString, description = ''))
									
						#output one sequence, for the current ST with all the metadata grouped
						if args.jgroup: phyloSeq.append(SeqRecord(Seq(stSeq,IUPAC.unambiguous_dna),id=bacterium+'_ST'+str(profileCode)+'_'+descriptionString, description = ''))
					else:
						for profileInstance in STmapper[profileCode]:
							metadataPointer = metadataJoinField if metadataJoinField in profileInstance else 'sampleID'
							
							phyloSeq.append(SeqRecord(Seq(stSeq,IUPAC.unambiguous_dna),id=bacterium+'_ST'+str(profileCode)+'_'+profileInstance[metadataPointer], description = ''))

				
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
								phyloSeq.append(SeqRecord(Seq(stSeq,IUPAC.unambiguous_dna),id=bacterium+'_ST'+str(profileCode)+'_'+str(prog)+'_'+descriptionString, description = ''))
					#output one sequence, for the current ST with all the metadata grouped
					if args.jgroup: phyloSeq.append(SeqRecord(Seq(stSeq,IUPAC.unambiguous_dna),id=bacterium+'_ST'+str(profileCode)+'_'+descriptionString, description = ''))
				else:
					for profileInstance in STmapper[profileCode]: 
						metadataPointer = metadataJoinField if metadataJoinField in profileInstance else 'sampleID'
						phyloSeq.append(SeqRecord(Seq(stSeq,IUPAC.unambiguous_dna),id=bacterium+'_ST'+str(profileCode)+'_'+profileInstance[metadataPointer], description = ''))
			
			SeqIO.write(phyloSeq,args.folder+'/merged/'+bacterium+'_sequences.fna', "fasta")

print "Color Legend:\n"+"-"*80
print "Alleles:"+'\t'+"[Known]"+'\t'+bcolors.OKBLUE+"[NEW]"+bcolors.ENDC+'\t'+bcolors.OKGREEN+"[NEW-RECURRING]"+bcolors.ENDC
print "Profiles:"+'\t'+bcolors.FAIL+"[Known]"+bcolors.ENDC+'\t'+bcolors.WARNING+"[NEW]"+bcolors.ENDC+'\t'+bcolors.OKGREEN+"[NEW*]"+bcolors.ENDC
print "New* profiles are composed by Known and Recurring alleles only"
print '-'*80


print "Completed! Have a nice day."
