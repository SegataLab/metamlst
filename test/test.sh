#Create Database with the sequences from MLST_sepidermidis.fasta"

../metaMLST-index.py -d sepidermidis.db -s MLST_sepidermidis.fasta

#Create Database with the typings from MLST_sepidermidis_types.txt"
../metaMLST-index.py -d sepidermidis.db -t MLST_sepidermidis_types.txt

#Generate a Bowtie2 index
../metaMLST-index.py -d sepidermidis.db -i bowtie_sepidermidis

#Map the fastq with Bowtie
bowtie2 --threads 4 --very-sensitive-local -a --no-unal -x bowtie_sepidermidis -U SRS013261_epidermidis.fastq | samtools view -bS - > SRS013261_epidermidis.bam;

#Run MetaMLST on a single sample
../metaMLST.py -d sepidermidis.db SRS013261_epidermidis.bam -o ./out/

#Type the STs
../metaMLST-merge.py -d sepidermidis.db ./out

