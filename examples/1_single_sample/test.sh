echo "
-------------------------------------------------------------------------------
          METAMLST TEST SCRIPT # 1 - ONE SAMPLE, DEFAULT DATABASE
-------------------------------------------------------------------------------

This test script executes MetaMLST on a single sample using the pre-made
database
-------------------------------------------------------------------------------

Input     : A small subset of HMP sample SRS013261, FASTQ format

Database  : Pre-made database: ../../metamlstDB_2017.db.

Detection : SRS013261: S. epidermidis ST-19

Output    : - Report file:
              - ./out/merged/sepidermidis_report.txt
            - Updated Typing list (STs):
              - ./out/merged/sepidermidis_ST.txt

-------------------------------------------------------------------------------"

#Generate a Bowtie2 index
echo "
-------------------------------------------------------------------------------
Executing: metaMLST-index.py -i bowtie_MmetaMLST ../../metamlstDB_2017.db
-------------------------------------------------------------------------------
"

../../metaMLST-index.py -i bowtie_MmetaMLST ../../metamlstDB_2017.db

#Map the fastq with Bowtie
echo "
-------------------------------------------------------------------------------
Executing: bowtie2 --threads 4 --very-sensitive-local -a --no-unal -x bowtie_MmetaMLST -U SRS013261_epidermidis.fastq | samtools view -bS - > SRS013261_epidermidis.bam;
-------------------------------------------------------------------------------
"

bowtie2 --threads 4 --very-sensitive-local -a --no-unal -x bowtie_MmetaMLST -U SRS013261_epidermidis.fastq | samtools view -bS - > SRS013261_epidermidis.bam;

#Run MetaMLST on a single sample
echo "
-------------------------------------------------------------------------------
Executing: metaMLST.py -d sepidermidis.db SRS013261_epidermidis.bam -o ./out/
-------------------------------------------------------------------------------
"

../../metaMLST.py -d ../../metamlstDB_2017.db SRS013261_epidermidis.bam -o ./out/

#Type the STs
echo "
-------------------------------------------------------------------------------
Executing: metaMLST-merge.py -d sepidermidis.db ./out
-------------------------------------------------------------------------------
"

../../metaMLST-merge.py -d ../../metamlstDB_2017.db ./out

