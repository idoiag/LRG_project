# LRG_project
Contains documentation about the LRG_project

Aim: LRG_ext extracts information about the introns and exons from a specified LRG identifier.
It also provides information about the GRCh37 and 38.

Structure:
The program runs several tests to test the input and if an output can be provided

  a) Tests (Control BRCA1)
    1. LRG_ID: check if the given gene has an allocated LRG ID
    2. Strain: check the direction of the strand. 
    3. Build_number: check the number of build provided. If more than one, check that coordinates are diff. 
    4. Build_coor: check that the start and end of the coordinate for each build are different
    5. Data_consistency: check that coordinates extracted are numbers
    6. Exon_coor: check that exons do not overlap
    7. Exon_number: check that we get the right number of exons
  
  b) modules
    1. Input name of gene 
    2. Extraction of LRG_identifier (if any)
    3. 
    2.
    3.
  
  
