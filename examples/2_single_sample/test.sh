echo "---- Executing: metaMLST-index.py -i bowtie_MmetaMLST ../../metamlstDB_2015.db"
#Create Database with the sequences from MLST_sepidermidis.fasta"
../../metaMLST-index.py -s MLST_sepidermidis.fasta sepidermidis.db 

echo "---- Executing: metaMLST-index.py -t MLST_sepidermidis_types.txt sepidermidis.db"
#Create Database with the typings from MLST_sepidermidis_types.txt"
../../metaMLST-index.py -t MLST_sepidermidis_types.txt sepidermidis.db

echo "---- Executing: metaMLST-index.py -i bowtie_sepidermidis sepidermidis.db"
#Generate a Bowtie2 index
../../metaMLST-index.py -i bowtie_sepidermidis sepidermidis.db

echo "---- Executing: bowtie2 --threads 4 --very-sensitive-local -a --no-unal -x bowtie_sepidermidis -U SRS013261_epidermidis.fastq | samtools view -bS - > SRS013261_epidermidis.bam;"
#Map the fastq with Bowtie
bowtie2 --threads 4 --very-sensitive-local -a --no-unal -x bowtie_sepidermidis -U SRS013261_epidermidis.fastq | samtools view -bS - > SRS013261_epidermidis.bam;

echo "---- Executing: metaMLST.py -d sepidermidis.db SRS013261_epidermidis.bam -o ./out/"
#Run MetaMLST on a single sample
../../metaMLST.py -d sepidermidis.db SRS013261_epidermidis.bam -o ./out/

echo "---- Executing: metaMLST-merge.py -d sepidermidis.db ./out"
#Type the STs
../../metaMLST-merge.py -d sepidermidis.db ./out

