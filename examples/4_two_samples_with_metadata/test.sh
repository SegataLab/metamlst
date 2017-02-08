echo "
-------------------------------------------------------------------------------
          METAMLST TEST SCRIPT # 4 - TWO SAMPLES, WITH METADATA
-------------------------------------------------------------------------------

This test script executes MetaMLST on two sample, using the pre-made
database, and detects the two species present.
-------------------------------------------------------------------------------

Input    : A small subset of HMP sample SRS013261, FASTQ format
           A small subset of HMP sample SRS015937, FASTQ format
           A metadata-file: test_metadata.txt, tab-separated format

Database : Pre-made database: ../../metamlstDB_2017.db


Detection : SRS013261: S. epidermidis ST-19
                       P. acnes ST-4

            SRS015937: S. epidermidis ST-100001 (New ST)

Output    : - Report files:
              - ./out/merged/sepidermidis_report.txt
              - ./out/merged/pacnes_report.txt

            - Updated Typing lists (STs):
              - ./out/merged/sepidermidis_ST.txt
              - ./out/merged/pacnes_ST.txt

-------------------------------------------------------------------------------"

#Generate a Bowtie2 index
echo "
-------------------------------------------------------------------------------
Executing: metamlst-index.py -i bowtie_MmetaMLST ../../metamlstDB_2017.db
-------------------------------------------------------------------------------
"
../../metamlst-index.py -i bowtie_MmetaMLST ../../metamlstDB_2017.db

#Map the fastq with Bowtie
echo "
-------------------------------------------------------------------------------
Executing: bowtie2 --threads 4 --very-sensitive-local -a --no-unal -x bowtie_MmetaMLST -U SRS015937_epidermidis.fastq | samtools view -bS - > SRS015937_epidermidis.bam
-------------------------------------------------------------------------------
"

bowtie2 --threads 4 --very-sensitive-local -a --no-unal -x bowtie_MmetaMLST -U SRS015937_epidermidis.fastq | samtools view -bS - > SRS015937_epidermidis.bam;

echo "
-------------------------------------------------------------------------------
Executing: bowtie2 --threads 4 --very-sensitive-local -a --no-unal -x bowtie_MmetaMLST -U SRS013261_epidermidis.fastq | samtools view -bS - > SRS013261_epidermidis.bam
-------------------------------------------------------------------------------
"

bowtie2 --threads 4 --very-sensitive-local -a --no-unal -x bowtie_MmetaMLST -U SRS013261_epidermidis.fastq | samtools view -bS - > SRS013261_epidermidis.bam;

#Run MetaMLST on a single sample
echo "
-------------------------------------------------------------------------------
Executing: metamlst.py -d sepidermidis.db SRS015937_epidermidis.bam -o ./out/
-------------------------------------------------------------------------------
"
../../metamlst.py -d ../../metamlstDB_2017.db SRS015937_epidermidis.bam -o ./out/

echo "
-------------------------------------------------------------------------------
Executing: metamlst.py -d sepidermidis.db SRS013261_epidermidis.bam -o ./out/
-------------------------------------------------------------------------------
"

../../metamlst.py -d ../../metamlstDB_2017.db SRS013261_epidermidis.bam -o ./out/

#Type the STs
echo "
-------------------------------------------------------------------------------
Executing: metamlst-merge.py -d sepidermidis.db --meta test_metadata.txt ./out
-------------------------------------------------------------------------------
"
../../metamlst-merge.py -d ../../metamlstDB_2017.db --meta test_metadata.txt ./out

