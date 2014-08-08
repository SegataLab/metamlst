import sys,os,subprocess,sqlite3,argparse
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


	

parser = argparse.ArgumentParser()
parser.add_argument("folder", help="folder containing .txt sametagenome files")
parser.add_argument("-d","--database", help="database", default="bsb.db")
parser.add_argument("--meta", help="metadata file (CSV)")
parser.add_argument("--idField", help="field containig the 'sampleID' value", default=0, type=int)

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
		if result: recs.append(str(result['recID']))

	return [(row['profileCode'],int((float(row['T'])/float(len(recs)))*100)) for row in e.execute("SELECT profileCode, COUNT(*) as T FROM profiles WHERE alleleCode IN ("+','.join(recs)+") GROUP BY profileCode HAVING T = (SELECT COUNT(*)  FROM profiles WHERE alleleCode IN ("+','.join(recs)+") GROUP BY profileCode ORDER BY COUNT(*) DESC LIMIT 1) ORDER BY T DESC")] if result else [(0,0)]
	
def sequenceExists(bacterium,sequence):
	
	e = conn.cursor()
	e.execute("SELECT 1 FROM alleles WHERE sequence = ? AND bacterium = ?",(str(sequence),bacterium))
	return (len(e.fetchall()) > 0)

def sequenceLocate(bacterium,sequence):
	
	e = conn.cursor()
	e.execute("SELECT alleleVariant FROM alleles WHERE sequence = ? AND bacterium = ?",(str(sequence),bacterium))
	return str(e.fetchone()['alleleVariant'])

cel = {}
 
for file in os.listdir(args.folder):
	
	if file.split('.')[-1] != 'txt': continue 
 
	
	for line in open(args.folder+'/'+file,'r'):
		organism = line.split()[0] 
		sampleName = line.split()[1] 
		genes =line.split()[2::] 
		
		if organism not in cel: cel[organism] = []
		cel[organism].append((dict((x.split('::')[0],(x.split('::')[1].upper(),x.split('::')[2])) for x in genes),sampleName))

for bacterium,bactRecord in cel.items():
	#bactRecord = ( {gene_allele: (sequence,accuracy) , ...} , sampleName ) 
	
	print '\r\n----------------------------------------------------------------------\r\n'+bacterium  
	 
	profil = open(args.folder+'/ST_'+bacterium+'.txt','w')
	
 	
	profil.write(bacterium+'\r\n')
	
	stringBase = {}
	
	oldProfiles = {}
	genesBase={}
	profilesBase={}
	encounteredProfiles={}
	isolates=[]
	newSequences = []
	
	cursor.execute("SELECT profileCode FROM profiles WHERE bacterium = ? ORDER BY profileCode DESC LIMIT 1",(bacterium,))

	lastProfile = int(cursor.fetchone()['profileCode'])
	lastGenes = dict((row['gene'],row['maxGene']) for row in cursor.execute("SELECT gene, MAX(alleleVariant) as maxGene FROM alleles WHERE bacterium = ? GROUP BY gene",(bacterium,)))
	
	#print "LAST PROFILE: ",lastProfile,"LAST GENES: ", lastGenes.items()
	
	#Get all "DB" profiles
	
	for row in cursor.execute("SELECT profileCode,gene,alleleVariant FROM profiles,alleles WHERE alleleCode = alleles.recID AND alleles.bacterium = ?",(bacterium,)):
		if row['profileCode'] not in oldProfiles: oldProfiles[row['profileCode']] = [0,{}]
		oldProfiles[row['profileCode']][1][row['gene']] = row['alleleVariant']
	
	
	
	profil.write('ST\t'+'\t'.join([x for x in sorted(lastGenes.keys())] )+'\r\n')
	print 'ST\t'+'\t'.join([x for x in sorted(lastGenes.keys())])+'\tHits'
	
	
	for bacteriumLine,sampleRecord in bactRecord:
		 
		#bacteriumLine is a dict. Key: genes, Values: sequences. No sequence = no news.
		#>> {gene_allele: (sequence,accuracy) , ...} 
		profileLine = {}
		newAlleles = []
		sum_of_accuracies = 0.0
		
		for geneLabel,(geneSeq,geneAccur) in bacteriumLine.items():
			
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
				
			elif geneSeq not in genesBase:
					## WE HAVE A NEW SEQUENCE
				geneNewAlleleNumber = str(lastGenes[geneName]+1)
				lastGenes[geneName]+=1
				geneNewLabel = geneOrganism+'_'+geneName+'_'+geneNewAlleleNumber
				genesBase[geneSeq] = geneNewLabel 
				#print "\t New Sequence for "+geneName+" -> "+geneNewLabel
				profileLine[geneName] = (geneNewAlleleNumber,1)
				newAlleles.append(geneName)
				newSequences.append(SeqRecord(Seq(geneSeq,IUPAC.unambiguous_dna),id=geneNewLabel, description = ''))
 
		
		meanAccuracy = sum_of_accuracies / float(len(bacteriumLine))
		if len(newAlleles) == 0:
		
		## Existent MLST profile -> Look for it (--> ISOLATES)
		
			#Tries to define an existing MLST profile with the alleles
			
			tryDefine = defineProfile([bacterium+'_'+k+'_'+v[0] for k,v in profileLine.items()])
			if tryDefine and tryDefine[0][1] == 100: 
			
				#encounteredProfiles[tryDefine[0][0]] = [profileLine,1,0]
				# TODO: to isolates and report as "non-new"
				print "\t | Existing MLST Profile", tryDefine[0][0]
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
			encounteredProfiles[lastProfile] = [profileLine,1,1]
			isolates.append((lastProfile,meanAccuracy,sampleRecord)) 
		
		
			
		
		# PROFILE LINE: dictionary['gene']: allele, hits, color
		#	-- ExistingCode | Effect --
		#           1       |  GREEN  -> New, with new alleles
		#           2       |  YELLOW -> New, with old allelese
		#                   |
		#                   | 
	#Old profiles
	
	SeqIO.write(newSequences,args.folder+'/sequences_'+bacterium+'.faa', "fasta")
	
	for profileCode,profileTrack in oldProfiles.items():
		profil.write(str(profileCode)+'\t'+'\t'.join([str(v) for k,v in sorted(profileTrack[1].items())])+'\r\n')
		
	#New profiles
		
	for profileID,(profile,hits,isNewProfile) in encounteredProfiles.items():
		
		if isNewProfile == 1: profileNumber = bcolors.OKGREEN+str(profileID)+bcolors.ENDC
		elif isNewProfile == 2: profileNumber = bcolors.WARNING+str(profileID)+bcolors.ENDC
		
		print profileNumber + '\t' + '\t'.join([bcolors.OKBLUE+str(v[0])+bcolors.ENDC if v[1] == 1 else bcolors.OKGREEN+str(v[0])+bcolors.ENDC if v[1] == 2 else str(v[0]) for k,v in sorted(profile.items())]) + '\t' + str(hits)
			
		profil.write(str(profileID)+'\t'+'\t'.join([str(v[0]) for k,v in sorted(profile.items())])+'\n')
	profil.close()
	
	if args.meta:
		isolafil = open(args.folder+'/isolate_'+bacterium+'.txt','w') 
		identifiers = {}
		p1line = 0
		for line in open(args.meta):
			if not p1line: 
				p1line=1
				keys = [str(x).strip() for x in line.split('\t')]
			else:
				l = line.strip().split('\t')
				identifiers[l[args.idField]] = dict((keys[i],l[i]) for i in range(0,len(keys)))
 
		 
		
		isolafil.write('ST\tAccuracy\t'+'\t'.join(keys)+'\n')
		
		for profileST,meanAccur,sampleName in isolates:
			if sampleName in identifiers:
				strl=[]
				for ky in keys:
					strl.append(identifiers[sampleName][ky])
				isolafil.write(str(profileST)+'\t'+str(round(meanAccur,2))+'\t'+'\t'.join(strl)+'\n')
			else:
				print "SAMPLE ",sampleName,"Not in identifiers"
				isolafil.write(str(profileST)+'\t'+str(round(meanAccur,2))+'\t'+str(sampleName)+'\n')
				
		isolafil.close()
		 
	
print "Completed"