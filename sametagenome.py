import argparse, re, os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import json
import math


parser = argparse.ArgumentParser()
parser.add_argument("bowfile", help="BOWTIE2 file containing the sequences")
#parser.add_argument("-i,--text", help="Information on the sample" action="store_true")
parser.add_argument("-o","--organism", help="target one species")
parser.add_argument("-p","--penalty", help="penalty for not found allele", default=100, type=int)
parser.add_argument("--alleles", help="number of top-fitting alleles to be shown", default=1000, type=int)
parser.add_argument("--minscore", help="minimum score to match (absolute val!)", default=50, type=int)
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
			else: ignoredReads = ignoredReads+1

for speciesKey,species in cel.items():
	for geneKey, geneInfo in species.items(): #geni
		maxLen = 0
		for geneInfoKey,geneValues in geneInfo.items(): #alleli passata 1
			if len(geneValues) > maxLen: maxLen = len(geneValues)

		for geneInfoKey,geneValues in geneInfo.items(): #alleli passata 2
			geneLen = len(geneValues)
			if geneLen == maxLen: cel[speciesKey][geneKey][geneInfoKey] = (sum(item for item in geneValues),maxLen)
			else: cel[speciesKey][geneKey][geneInfoKey] = (sum(item for item in geneValues) + (maxLen-geneLen)*args.penalty ,geneLen)
		
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

for speciesKey,species in cel.items():
	dfil.write(speciesKey)
	dfil.write('\r\n{')
	for geneKey, geneInfo in species.items(): #geni
		dfil.write('\t'+geneKey+'\t'+repr((geneMaxiMin[speciesKey+geneKey])[::-1])+'\r\n')
		dfil.write('\t{\r\n')
		for geneInfoKey,geneValues in sorted(geneInfo.items(), key= lambda x: x[1]): #alleli
			dfil.write('\t\t'+geneKey + geneInfoKey+'\t\t'+repr(geneValues[0])+'\t\t'+repr(geneValues[1])+'\r\n')
		dfil.write('\t}\r\n')
	dfil.write('}\r\n')
dfil.close()	
###################################### OUT FILE 2

dfil = open(fileName+'/'+fileName+'.out2','w') 
dfil.write("# intest #")

for speciesKey,species in cel.items():
	#dfil.write(speciesKey)
	#dfil.write('\r\n{')
	for geneKey, geneInfo in species.items(): #geni
		#dfil.write('\t'+geneKey+'\t'+repr((geneMaxiMin[speciesKey+geneKey])[::-1])+'\r\n')
		#dfil.write('\t{\r\n')
		max = float('inf')
		
		minv = min([x[0] for x in geneInfo.values()])
		
		geneInfo = dict([(k,v) for k,v in geneInfo.items() if v[0] == minv])
		for geneInfoKey,geneValues in sorted(geneInfo.items(), key= lambda x: x[1]): #allele top
			dfil.write(str(speciesKey)+'\t'+str(geneKey)+'\t'+str(geneInfoKey)+'\t'+str(geneValues[0])+'\t'+str(geneValues[1])+'\t'+str(geneValues[0]/geneValues[1])+'\t\r\n')
		#dfil.write('\t}\r\n')
	#dfil.write('}\r\n')
dfil.close()	

for sequenceKey,sequenceRecord in sequenceBank.items():
	dfil = open(fileName+'/'+sequenceKey+'.fastq','a')
	for sequenceSpec,(sequence,quality) in sequenceRecord.items():
		dfil.write('@'+sequenceSpec+'\r\n')  
		dfil.write(sequence+'\r\n')  
		dfil.write('+\r\n')  
		dfil.write(quality+'\r\n')  
	dfil.close()

print fileName+' Completed'
#print cel
