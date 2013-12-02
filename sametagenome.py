import argparse, re, os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import json, math,sqlite3


parser = argparse.ArgumentParser()
parser.add_argument("bowfile", help="BOWTIE2 file containing the sequences")
#parser.add_argument("-i,--text", help="Information on the sample" action="store_true")
parser.add_argument("-o","--organism", help="target one species")
parser.add_argument("-p","--penalty", help="penalty for not found allele", default=100, type=int)
parser.add_argument("--alleles", help="number of top-fitting alleles to be shown", default=1000, type=int)
parser.add_argument("--minscore", help="minimum score to match (absolute val!)", default=50, type=int)
parser.add_argument("--mincoverage", help="minimum coverage to produce FASTQ assembler files", default=5, type=int)
parser.add_argument("--min_read_len", help="minimum length for bowtie2 reads", default=90, type=int)
parser.add_argument("--graphic", help="generate graphics", action="store_true")
args=parser.parse_args()

fil = open(args.bowfile,'r')
cel = {}
geneMaxiMin = {}
sequenceBank = {}
ignoredReads = 0

for line in fil:
	if(line[0] != '@'): 
		read = line.split('\t')
		readCode = read[0]
		species = re.findall('[a-zA-Z]*_',read[2])[0].replace('_','')
		gene = re.findall('_[a-zA-Z_]*',read[2])[0].replace('_','')
		allele = re.findall('[0-9]*$',read[2])[0]
		score = abs((int)((read[11]).split(':')[2]))
		sequence = read[9]
		quality = read[10]
		#quality = [(ord(x)-ord("!")) for x in read[10]] ???
		
		if (args.organism and species == args.organism) or not args.organism:
			if score <= args.minscore and len(sequence) >= args.min_read_len:
				if species not in cel: cel[species] = {} #new species
			
				if gene not in cel[species]: #new gene
					cel[species][gene] = {}
					geneMaxiMin[species+gene] = [0,float("inf")]
					sequenceBank[species+'_'+gene] = {}
					
				if allele not in cel[species][gene]: #new allele
					cel[species][gene][allele] = []
				
				cel[species][gene][allele].append(score)
				#sequenceBank[species+'_'+gene].append(SeqRecord(Seq(sequence, IUPAC.unambiguous_dna), id = readCode, description = ''))
				#sequenceBank[species+'_'+gene].append((readCode,sequence,quality))
				sequenceBank[species+'_'+gene][readCode]=(sequence,quality)
				#print readCode
			else: ignoredReads = ignoredReads+1
			
for speciesKey,species in cel.items():
	for geneKey, geneInfo in species.items(): #geni
	
		maxLen = max([len(x) for x in geneInfo.values()])

		for geneInfoKey,geneValues in geneInfo.items(): #alleli passata 2
			geneLen = len(geneValues)
			localScore = sum(item for item in geneValues)
			
			averageScore = float(localScore)/float(geneLen)
			
			if geneLen != maxLen:
				localScore = localScore + (maxLen-geneLen)*args.penalty
			
			cel[speciesKey][geneKey][geneInfoKey] = (localScore,maxLen,round(averageScore,1)) 
		
			if cel[speciesKey][geneKey][geneInfoKey][0] > geneMaxiMin[speciesKey+geneKey][0]:
				geneMaxiMin[speciesKey+geneKey][0] = cel[speciesKey][geneKey][geneInfoKey][0]
				
			if cel[speciesKey][geneKey][geneInfoKey][0] < geneMaxiMin[speciesKey+geneKey][1]:
				geneMaxiMin[speciesKey+geneKey][1] = cel[speciesKey][geneKey][geneInfoKey][0]

fileName = (args.bowfile.split('/'))[-1].split('.')[0]
os.mkdir(fileName)
				
if args.graphic:
	jsonO = {'nodes': [],'links': []}
	tracker = {'genes':{},'alleles':{},'species':{},'speciesRight' : {}}
	idP = 0;

	for species, speciesInfo in cel.items():
		jsonO['nodes'].append({'name':species})
		tracker['species'][species] = idP
		idP = idP+1
		#jsonO['nodes'].append({'name':species})
		#tracker['speciesRight'][species] = idP
		#idP = idP+1
		for gene, geneInfo in speciesInfo.items(): #geni della Specie
			jsonO['nodes'].append({'name':gene})
			tracker['genes'][gene] = idP
			idP = idP+1
			valSum = 0;
			i=0
			#for allele,alleleRecord in geneInfo.items():
			#for allele,alleleRecord in sorted(geneInfo.items(), key = lambda x: x[1]):
			for allele,alleleRecord in sorted(geneInfo.items(), key= lambda x: x[1]): #alleli
				if i < args.alleles:
					jsonO['nodes'].append({'name':gene+allele})
					tracker['alleles'][allele] = idP
					idP = idP+1
					range = (geneMaxiMin[species+gene][0] - geneMaxiMin[species+gene][1])
					if range == 0 : range = 1
					
					#val = geneMaxiMin[species+gene][0] - (alleleRecord[0]-geneMaxiMin[species+gene][1])
					val = round((1-(float)(alleleRecord[0]-geneMaxiMin[species+gene][1]) / (float)(range)),4)*100
					
					#print gene,allele,'\t\t',alleleRecord[0],geneMaxiMin[species+gene],range,val
					jsonO['links'].append({"source":tracker['genes'][gene],"target":tracker['alleles'][allele],"value":round(val)})
					#jsonO['links'].append({"source":tracker['alleles'][allele],"target":tracker['speciesRight'][species],"value":round(val)})
					valSum = valSum + val
				i=i+1
				
			jsonO['links'].append({"source":tracker['species'][species],"target":tracker['genes'][gene],"value":round(valSum)})
			#fine allele
		#fine gene
	#fine specie
	
	#stampa
	tfil = open('template.txt','r')
	dfil = open(fileName+'/'+fileName+'.htm','w')
	dfil.write('<script>var dataSet = '+json.dumps(jsonO)+';</script>'+tfil.read())
	dfil.close()
	tfil.close()

dfil = open(fileName+'/'+fileName+'.out','w')
dfil.write("SAMPLE:\t\t\t\t\t"+args.bowfile+'\r\n')
dfil.write("PENALTY:\t\t\t\t"+repr(args.penalty)+'\r\n')
dfil.write("MAX-ALLELES:\t\t\t"+repr(args.alleles)+'\r\n')
dfil.write("MIN-THRESHOLD SCORE:\t"+repr(args.minscore)+'\r\n')
dfil.write("IGNORED:\t\t\t\t"+repr(ignoredReads)+' BAM READS\r\n\r\n------------------------------  RESULTS ------------------------------\r\n')

# cel ['haemophilus']['haemophilus_gene'][allele] = (a, b, c) --> a: summed score, b = maxLen, c = a/b

###################################### OUT FILE 1 (LOG)
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
###################################### OUT FILE 2 (USEFUL)


dfil = open(fileName+'/'+fileName+'.out2','w') 
dfil.write("# intest #\r\n\r\n")

print "Sample Contains:\r\n-----------------------------------"

for speciesKey,species in cel.items():
	#dfil.write(speciesKey)
	#dfil.write('\r\n{')
	
	conn = sqlite3.connect('bsb.db')
	conn.row_factory = sqlite3.Row
	c = conn.cursor()
	d = conn.cursor()
	tVar = dict([(row['geneName'],0) for row in  c.execute("SELECT geneName FROM genes WHERE bacterium = ?",(speciesKey,))])
	
	#GENE PRESENCE 
	for sk in species.keys():
		tVar[sk] = 1
	vals = sum([t for t in tVar.values()])
	#
	
	# NOT ENOUGH GENES
	if float(vals)/float(len(tVar)) < 0.8:
		dfil.write('#'+str(speciesKey)+'\t SKIPPED: Not enough genes ('+str(vals)+' / '+str(len(tVar))+')\r\n')
		print "\033[91m", speciesKey,"\033[0m\t: ", str(vals)+' genes out of '+str(len(tVar))+' MLST targets'
		print "\t\t   Missing: \033[91m"+', '.join([sk for (sk,v) in tVar.items() if v == 0])+'\033[0m'
		
	# ENOUGH GENES
	else:
		print "\033[92m", speciesKey,"\033[0m\t: ", str(vals)+' genes out of '+str(len(tVar))+' MLST targets'
		print "\t\t   Missing: \033[91m"+', '.join([sk for (sk,v) in tVar.items() if v == 0])+'\033[0m'
		
		print "\r\n  "+"Gene".ljust(7)+"Coverage".rjust(10)+"Score".rjust(6)+"Hits".rjust(5)+" Allele(s)".ljust(40)
		
		profileTrack = []
		for geneKey, geneInfo in species.items(): #geni
		
			minValue = min([avg for (val,leng,avg) in geneInfo.values()])
			
			aElements = {}
			tmp=[]
			for k,(val,leng,avg) in geneInfo.items():
				if avg == minValue:
					aElements[k]=(val,leng,avg)
					tmp.append(k)
			tmp = ",".join(tmp)
			
			profileTrack = profileTrack + [str(row['recID']) for row in c.execute("SELECT recID FROM alleles WHERE bacterium = ? AND gene = ? AND alleleVariant IN ("+tmp+")",(speciesKey,geneKey))]
			
			dfil.write(str(speciesKey)+'\t'+str(geneKey)+'\t'+str(minValue)+'\t'+repr([geneKey+k for k in aElements])+'\r\n')
			
			sequenceKey = speciesKey+'_'+geneKey
			c.execute("SELECT LENGTH(sequence) as L FROM alleles WHERE bacterium = ? AND gene = ? ORDER BY L DESC LIMIT 1", (speciesKey,geneKey))
			genL = c.fetchone()['L']
		
			coverage = sum([len(x) for (x,q) in sequenceBank[sequenceKey].values()])
			if coverage >= args.mincoverage * genL:
				fqfil = open(fileName+'/'+sequenceKey+'.fastq','a')
				color = "\033[92m"
				for sequenceSpec,(sequence,quality) in sequenceBank[sequenceKey].items():
					fqfil.write('@'+sequenceSpec+'\r\n')  
					fqfil.write(sequence+'\r\n')  
					fqfil.write('+\r\n')  
					fqfil.write(quality+'\r\n') 
				fqfil.close()
			else:  color = "\033[93m"
			
			print "  "+color+geneKey.ljust(7)+"\033[0m"+str(round(float(coverage)/float(genL),2)).rjust(10)+"\033[95m"+str(minValue).rjust(6)+str(aElements.itervalues().next()[1]).rjust(5)+"\033[0m"+"\033[94m",tmp.ljust(40)+"\033[0m"
			
		print ""
		
		for row in c.execute("SELECT profileCode, COUNT(*) as T FROM profiles WHERE alleleCode IN ("+','.join(profileTrack)+") GROUP BY profileCode HAVING T = (SELECT COUNT(*)  FROM profiles WHERE alleleCode IN ("+','.join(profileTrack)+") GROUP BY profileCode ORDER BY COUNT(*) DESC LIMIT 1) ORDER BY T DESC"):
			matchScore = str(round(float(row['T']) / float(len(tVar)),4)*100)+' %'
			print ("  MLST PROFILE "+str(row['profileCode'])).ljust(23)+str("MATCH: "+matchScore).rjust(10)
			
			v = [riw['sequence'] for riw in d.execute("SELECT sequence FROM profiles,alleles WHERE alleleCode = alleles.recID AND profiles.bacterium = ? AND profileCode = ?",(speciesKey,row['profileCode']))]
dfil.close()	
conn.close() 

print fileName+'\t\033[92mCompleted\033[0m'
#print cel
