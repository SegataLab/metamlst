echo "
-------------------------------------------------------------------------------
          METAMLST TEST SCRIPT # 2 - ONE SAMPLE, CUSTOM DATABASE
-------------------------------------------------------------------------------

This test script executes MetaMLST on a single sample using a custom database
generated from your sequences and typing tables.
-------------------------------------------------------------------------------

Input     : A small subset of HMP sample SRS013261, FASTQ format

Database  : A CUSTOM database specific for S. epidermidis made from:
               - MLST_sepidermidis.fasta     (MLST sequences)
               - MLST_sepidermidis_types.txt (tab-separated typing file)

Detection : SRS013261: S. epidermidis ST-19
          
Output    : - Report files:
              - ./out/merged/sepidermidis_report.txt

            - Updated Typing lists (STs):
              - ./out/merged/sepidermidis_ST.txt

-------------------------------------------------------------------------------"

echo "
-------------------------------------------------------------------------------
Executing: metamlst-index.py -s MLST_sepidermidis.fasta sepidermidis.db
-------------------------------------------------------------------------------
"

#Create Database with the sequences from MLST_sepidermidis.fasta"
../../metamlst-index.py -s MLST_sepidermidis.fasta sepidermidis.db 

echo "
-------------------------------------------------------------------------------
Executing: metamlst-index.py -t MLST_sepidermidis_types.txt sepidermidis.db
-------------------------------------------------------------------------------
"

#Create Database with the typings from MLST_sepidermidis_types.txt"
../../metamlst-index.py -t MLST_sepidermidis_types.txt sepidermidis.db

echo "
-------------------------------------------------------------------------------
Executing: metamlst-index.py -i bowtie_sepidermidis sepidermidis.db
-------------------------------------------------------------------------------
"

#Generate a Bowtie2 index
../../metamlst-index.py -i bowtie_sepidermidis sepidermidis.db

echo "
-------------------------------------------------------------------------------
Executing: bowtie2 --threads 4 --very-sensitive-local -a --no-unal -x bowtie_sepidermidis -U SRS013261_epidermidis.fastq | samtools view -bS - > SRS013261_epidermidis.bam;"
#Map the fastq with Bowtie
bowtie2 --threads 4 --very-sensitive-local -a --no-unal -x bowtie_sepidermidis -U SRS013261_epidermidis.fastq | samtools view -bS - > SRS013261_epidermidis.bam;

echo "
-------------------------------------------------------------------------------
Executing: metamlst.py -d sepidermidis.db SRS013261_epidermidis.bam -o ./out/
-------------------------------------------------------------------------------
"

#Run MetaMLST on a single sample
../../metamlst.py -d sepidermidis.db SRS013261_epidermidis.bam -o ./out/

echo "
-------------------------------------------------------------------------------
Executing: metamlst-merge.py -d sepidermidis.db ./out
-------------------------------------------------------------------------------
"

#Type the STs
../../metamlst-merge.py -d sepidermidis.db ./out

