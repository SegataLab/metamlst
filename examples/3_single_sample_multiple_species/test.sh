echo "
-------------------------------------------------------------------------------
          METAMLST TEST SCRIPT # 3 - ONE SAMPLE, MULTIPLE SPECIES
-------------------------------------------------------------------------------

This test script executes MetaMLST on a single sample, using the pre-made
database, and detects the two species present.
-------------------------------------------------------------------------------

Input    : A small subset of HMP sample SRS013261, FASTQ format

Database : Pre-made database: ../../metamlstDB_2017.db.

Detection : SRS013261: S. epidermidis ST-19
                       P. acnes ST-4
          
Output    : - Report files:
              - ./out/merged/sepidermidis_report.txt
              - ./out/merged/pacnes_report.txt

            - Updated Typing lists (STs):
              - ./out/merged/sepidermidis_ST.txt
              - ./out/merged/pacnes_ST.txt

-------------------------------------------------------------------------------"

echo "
-------------------------------------------------------------------------------
Executing: metamlst-index.py -i bowtie_sepidermidis ../../metamlstDB_2017.db
-------------------------------------------------------------------------------
"

#Generate a Bowtie2 index
../../metamlst-index.py -i bowtie_sepidermidis ../../metamlstDB_2017.db

echo "
-------------------------------------------------------------------------------
Executing: bowtie2 --threads 4 --very-sensitive-local -a --no-unal -x bowtie_sepidermidis -U SRS013261_epidermidis.fastq | samtools view -bS - > SRS013261_epidermidis.bam;
-------------------------------------------------------------------------------
"

#Map the fastq with Bowtie
bowtie2 --threads 4 --very-sensitive-local -a --no-unal -x bowtie_sepidermidis -U SRS013261_epidermidis.fastq | samtools view -bS - > SRS013261_epidermidis.bam;

echo "
-------------------------------------------------------------------------------
Executing: metamlst.py -d ../../metamlstDB_2017.db SRS013261_epidermidis.bam -o ./out/
-------------------------------------------------------------------------------
"

#Run MetaMLST on a single sample
../../metamlst.py -d ../../metamlstDB_2017.db SRS013261_epidermidis.bam -o ./out/

echo "
-------------------------------------------------------------------------------
Executing: metamlst-merge.py -d ../../metamlstDB_2017.db ./out
-------------------------------------------------------------------------------
"

#Type the STs
../../metamlst-merge.py -d ../../metamlstDB_2017.db ./out

