echo "
-------------------------------------------------------------------------------
          METAMLST TEST SCRIPT # 5 - MANY SAMPLES (Phylogenetic Analysis)
-------------------------------------------------------------------------------

This test script executes MetaMLST on 27 samples, using the pre-made
database. The script detects P. acnes and reports files suitable for 
Phylogenetic Analysis (trees and Minimum Spanning Trees)
-------------------------------------------------------------------------------

Input     : 27 files from the HMP (filtered reads matching to P. acnes)

Database  : Pre-made database: ../../metamlstDB_2017.db

Detection : 13 known STs, 14 new STs (ID >100000)

Output    : - Report file: ./out/merged/pacnes_report.txt
            - ST file: ./out/merged/pacnes_ST.txt 
            - Sequences file: ./out/merged/pacnes_sequences.txt

-------------------------------------------------------------------------------"

#Unzipping the archive
unzip 5_phylogenetic_analysis_files.zip

#Generate a Bowtie2 index
echo "
-------------------------------------------------------------------------------
Executing: metaMLST-index.py -i bowtie_MmetaMLST ../../metamlstDB_2017.db
-------------------------------------------------------------------------------
"
../../metaMLST-index.py -i bowtie_MmetaMLST ../../metamlstDB_2017.db


for i in ./*.fastq; do
   
  bs=$(basename $i)
  echo "
-------------------------------------------------------------------------------
Executing: bowtie2 --threads 4 --very-sensitive-local -a --no-unal -x bowtie_MmetaMLST -U "$i" | samtools view -bS - > "${bs%%.*}".bam;
-------------------------------------------------------------------------------
" 
  #Map the fastq with Bowtie
  bowtie2 --threads 4 --very-sensitive-local -a --no-unal -x bowtie_MmetaMLST -U $i | samtools view -bS - > ${bs%%.*}.bam;
  
  echo "
-------------------------------------------------------------------------------
Executing: metaMLST.py -d ../../metamlstDB_2017.db "${bs%%.*}".bam --filter pacnes -o ./out/
-------------------------------------------------------------------------------
" 
  #Run MetaMLST
  ../../metaMLST.py -d ../../metamlstDB_2017.db ${bs%%.*}.bam --filter pacnes -o ./out/
done


#Type the STs
echo "
-------------------------------------------------------------------------------
Executing: metaMLST-merge.py -d ../../metamlstDB_2017.db --meta test_metadata.txt ./out
-------------------------------------------------------------------------------
"
../../metaMLST-merge.py -d ../../metamlstDB_2017.db --meta sample_metadata.txt --outseqformat A ./out




