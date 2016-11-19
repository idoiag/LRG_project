/Text within slashes is temp and will be deleted in a final version/

##############  README ############

#### BACKGROUND /VF/
/Examples: 
Why are LRGs are important, use, challenge
How can an xml parser could help the diagnostic service,
Etc… /


#### USAGE/INSTALLATION






#### FEATURES /VF + IGP/

The program runs several tests to test the input and if an output can be provided
  
    a) Tests (Control BRCA1)
    1. LRG_ID: check if the given gene has an allocated LRG ID
     2. Strain: check the direction of the strand. - IP
     3. Build_number: check the number of build provided. If more than one, check that coordinates are diff. 
     4. Build_coor: check that the start and end of the coordinate for each build are different
     5. Data_consistency: check that coordinates extracted are numbers
     6. Exon_coor: check that exons do not overlap
     7. Exon_number: check that we get the right number of exons
     8. Check the schema version - DONE
     10. Build_number: check the 
     9. Check the LRG exists and is a readable file - DONE but it might be superfluous. There is already a CL in-built system
     
     
    b) modules
      1. Access the LRG-sequence.org with an specific Gene name 
      2. Extraction of information
      3. Define the genomic coordinates for ezch exos regarding GRCh37 and GRCh38
      4. Compare coordnates
 

#### STRUCTURE /IGP/

…


Note: Refer to Developers file for further information





#### SCOPE /VF/
This software has been tested with a sample of .xml files obtained from the LRG website. We have tried to include the main genes analysed in diagnostic labstories, such as:

/Here we include a list of the genes tested, the gene name and its main impact/

LRG_292 – BRCA1 – Decisive in the diagnosis of Breast cancer
LRG_ ….
/N.B. LRG.xml files included for testing purposes include: LRG_1 (Gene COL1A1) LRG_9 (Gene SDHD) Pending approval, not public yet LRG_214 (Gene NF1) Two transcripts LRG_292 (Gene BRCA1) Reverse strand/





#### VERSIONING /IGP/

- Extracting background information about the gene: v1
- Making the code modular: v2
- Adding examples of initial and final test: v2.1
- Extracting exon, transcript and protein coordinates: v3
- Output: it creates a csv and tab separate text file: v3.1
- Dealing with more than one transcript v.4
- Adding disclaimer v.4.1
- Creation of a bed file v.5
- Organisation of input and outputs into folders v.5.1
- Automatic creation of output file name. Code clean, comments added v.5.2
- Bug resolved for LRG_292 and LRG 62 (Prot coordinates for exon 1 were missing) v.5.3


#### DISCLAIMER /IGP/
Please cite this software as: 'Gomez-Paramio, I. and Fryer, V. (2016), 'lrgext', Software, Faculty of Medicine and Human Sciences, The University of Manchester.' or successor references as defined by the authors.

For further information refer to the copyright


     
 

