# LRG_project
Contains documentation about the LRG_project

Aim: LRG_ext extracts information from the LRG record of a gene. This includes: 

- LRG ID
- Gene
- Chromosome number 
- strand
- Exons number 
- Exon number +  coordinates + genomic coordinates reg. GRCh37 and GRCh38)
- Compare coordinates between builds

Structure:
The program runs several tests to test the input and if an output can be provided

  a) Tests (Control BRCA1)
    1. Check the schema version - DONE
    2. Strain: check the direction of the strand  -- IP
    3. Build_number: check the 
    
    
    er of exons
    8. Check the LRG exists and is a readable file - DONE but it might be superfluous. There is already a CL in-built system
   
  
  b) modules
    1. Access the LRG-sequence.org with an specific Gene name 
    3. Extraction of LRG_identifier (if any), Gene, Chromosome number, strand, exons numbers
    4. Define the genomic coordinates for each exons regarding GRCh37 and GRCh38
    5. Compare coordinates

N.B. LRG.xml files included for testing purposes include:
LRG_1 (Gene COL1A1)
LRG_9 (Gene SDHD) Pending approval, not public yet
LRG_214 (Gene NF1) Two transcripts
LRG_292 (Gene BRCA1) Reverse strand
  
