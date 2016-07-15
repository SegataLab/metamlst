#Create Database with the sequences from MLST_sepidermidis.fasta and the typings"

../metaMLST-index.py -s MLST_sepidermidis.fasta -t MLST_sepidermidis_types.txt sepidermidis.db 

#Generate a Bowtie2 index
../metaMLST-index.py -i bowtie_sepidermidis sepidermidis.db

#Map the fastq with Bowtie
bowtie2 --threads 4 --very-sensitive-local -a --no-unal -x bowtie_sepidermidis -U SRS015937_epidermidis.fastq | samtools view -bS - > SRS015937_epidermidis.bam;
bowtie2 --threads 4 --very-sensitive-local -a --no-unal -x bowtie_sepidermidis -U SRS013261_epidermidis.fastq | samtools view -bS - > SRS013261_epidermidis.bam;


#Run MetaMLST on a single sample
../metaMLST.py -d sepidermidis.db SRS015937_epidermidis.bam -o ./out/
../metaMLST.py -d sepidermidis.db SRS013261_epidermidis.bam -o ./out/

#Type the STs
../metaMLST-merge.py -d sepidermidis.db --meta test_metadata.txt ./out

