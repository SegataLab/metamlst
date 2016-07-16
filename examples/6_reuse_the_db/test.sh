echo "
-------------------------------------------------------------------------------
          METAMLST TEST SCRIPT # 6 - ONE SAMPLE WITH DATABASE UPDATE
-------------------------------------------------------------------------------

This test script executes MetaMLST on a sample harboring a NEW ST (ID > 100000)

Then, MetaMLST-index is run on the resulting updated STs table and sequences 
file to build a custom database including the new type detected (ST 100001).

The pipeline is then run again on the same sample, but with the custom database
and the sample's S. epidermidis, previously detected as new, is now marked as
known (i.e. previously encountered). 

-------------------------------------------------------------------------------

Input     : A small subset of HMP sample SRS015937, FASTQ format
            A metadata-file: test_metadata.txt, tab-separated format

Database  : Pre-made database: ../../metamlstDB_2015.db


Detection : SRS015937: S. epidermidis ST-100001 (New ST, first pass)
            SRS015937: S. epidermidis ST-100001 (Known ST, second pass)

Output    : - Report files (first pass):
              - ./out/merged/sepidermidis_report.txt

            - Updated Typing lists (STs):
              - ./out/merged/sepidermidis_ST.txt

            - Sequences file (contains the loci-sequences):
              - ./out/merged/sepidermidis_sequences.fna

            - Report files (second pass):
              - ./NEW_out/merged/sepidermidis_report.txt

            - Updated Typing lists (STs):
              - ./NEW_out/merged/sepidermidis_ST.txt
-------------------------------------------------------------------------------"

#Generate a Bowtie2 index
echo "
-------------------------------------------------------------------------------
Executing: metaMLST-index.py -i bowtie_MmetaMLST ../../metamlstDB_2015.db
-------------------------------------------------------------------------------
"
../../metaMLST-index.py -i bowtie_MmetaMLST ../../metamlstDB_2015.db

#Map the fastq with Bowtie
echo "
-------------------------------------------------------------------------------
Executing: bowtie2 --threads 4 --very-sensitive-local -a --no-unal -x bowtie_MmetaMLST -U SRS015937_epidermidis.fastq | samtools view -bS - > SRS015937_epidermidis.bam
-------------------------------------------------------------------------------
"

bowtie2 --threads 4 --very-sensitive-local -a --no-unal -x bowtie_MmetaMLST -U SRS015937_epidermidis.fastq | samtools view -bS - > SRS015937_epidermidis.bam;

#Run MetaMLST on a single sample
echo "
-------------------------------------------------------------------------------
Executing: metaMLST.py -d sepidermidis.db SRS015937_epidermidis.bam -o ./out/
-------------------------------------------------------------------------------
"
../../metaMLST.py -d ../../metamlstDB_2015.db SRS015937_epidermidis.bam -o ./out/

#Type the STs
echo "
-------------------------------------------------------------------------------
Executing: metaMLST-merge.py -d ../../metamlstDB_2015.db --outseqformat B+ ./out
-------------------------------------------------------------------------------
"
../../metaMLST-merge.py -d ../../metamlstDB_2015.db --outseqformat B+ ./out
 

echo "
-------------------------------------------------------------------------------
Executing: metaMLST-index.py-s out/sepidermidis_sequences.fna -t out/sepidermidis_ST.txt new_metamlstDB.db
-------------------------------------------------------------------------------
"
echo "#sepidermidis|Staphylococcus epidermidis" > out/merged/update_sepidermidis_ST.txt
cat out/merged/sepidermidis_ST.txt >> out/merged/update_sepidermidis_ST.txt
../../metaMLST-index.py -s out/merged/sepidermidis_sequences.fna -t out/merged/update_sepidermidis_ST.txt new_metamlstDB.db

echo "
-------------------------------------------------------------------------------
Executing: metaMLST-index.py -i new_bowtie_MmetaMLST new_metamlstDB.db
-------------------------------------------------------------------------------
"
../../metaMLST-index.py -i new_bowtie_MmetaMLST new_metamlstDB.db

echo "
-------------------------------------------------------------------------------
Executing: bowtie2 --threads 4 --very-sensitive-local -a --no-unal -x new_bowtie_MmetaMLST -U SRS015937_epidermidis.fastq | samtools view -bS - > SRS015937_epidermidis.bam;
-------------------------------------------------------------------------------
"
bowtie2 --threads 4 --very-sensitive-local -a --no-unal -x new_bowtie_MmetaMLST -U SRS015937_epidermidis.fastq | samtools view -bS - > SRS015937_epidermidis.bam;

echo "
-------------------------------------------------------------------------------
Executing: metaMLST.py -d new_metamlstDB.db SRS015937_epidermidis.bam -o ./NEW_out/
-------------------------------------------------------------------------------
"
../../metaMLST.py -d new_metamlstDB.db SRS015937_epidermidis.bam -o ./NEW_out/

echo "
-------------------------------------------------------------------------------
Executing: metaMLST-merge.py -d new_metamlstDB.db ./NEW_out
-------------------------------------------------------------------------------
"
../../metaMLST-merge.py -d new_metamlstDB.db ./NEW_out