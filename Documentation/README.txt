/Text within slashes is temp and will be deleted in a final version. IP = In progress. Done = Only minor amendments/addition/updates expected/

##############  README ############

#### BACKGROUND /VF, IP/
/Examples: 
VF - Why LRGs are important, use, challenges
VF - How can an xml parser help the diagnostic service
LRG records are available online at www.lrg-sequence-org. These records contains a lot of data about each LRG and it is useful for the bioinformatician to be able to access these records and extract a subset of data that they and other clinical scientists within their department may find useful. Developing an xml parser that deals specifically with LRG records means that the exact data need can be extracted very quickly and consistantly each time the request is made.
Many laboratories currently store the LRG xml files on a local folder (which will decrease the processing time of the script slightly as it won't have to send requests via the internet). This method is also more secure as all the data processing can occur behind the NHS firewall. The disadvantage of this system is that the stored .xml files may not be of the most up-to-date information on the website, but this can be minimised by implemeting a regular routine for pulling .xml files from the site. The script itself compiles a list of all LRG xml files in the LRGs directory at the start of the script each time it is run therefore the script will be using the most up-to-date records that have been saved within the LRGs directory. This means that if additional or amended xml files are added to the directory, the script will use these and it won't fall over because it doesn't recognise a new LRG name or if one is missing that had previously been processed.
VF - How can bioinformaticians help > great idea to have a parser
Etc. /


#### USAGE/INSTALLATION /IGP, IP, pending on features being added/

lrgext does not require installation and can be run from the command line by providing the gene name as an argument (default extension is .xml). This script was originally using the LRG name as the argument, but as most users of the script will be most likely to know the gene name before the LRG, this script can be used to find out a) if an LRG exists (and has been saved within the LRGs folder) and b) what the LRG id is for a given gene name.
lrgext will extract build, gene, transcript, exon information from a file in LRG format 
producing different output files:
    - cvs file: provides the exon, transcript and protein coordinates separated by comma /IGP/
    - tab separated txt file: provides the exon, transcrip and protein coordinates separated by tab /IGP/
    - bed file: provides the local coordinates and their equivalent genomic coordinates for each exon /IGP/
    /- build mapper: provides the differences between the builds, highlighting if variants lie in intronic/exonic regions /VF//
	
/We  need to produce another doc highlighting the differences between the builds (VF) /

#### FEATURES /IGP + VF/

/Background:VF 
Example: sell our program highlighting that it is robust, efficient and its code could be reused since it was built in modules 
(see features below).
Also, it would be good if we could provide examples of having it built following the principles of defensive programming 
(i.e. using tests, see testing secton) but it will be important to build a test using assert. Ideally we could include a section of how we deal with error handling,


1.Testing: 
Testing is important because bla, bla... / 2-3 lines, VF/

/This is a list of tests we initally suggested. We also proposed to use BRCA as a control. 1-3 lines to describe each
function, IGP  + VF as required/

	1. LRG_ID: check if the given gene has an allocated LRG ID. /VF, Done
	2. Check the LRG exists and is a readable file. /IGP, DONE /
	3. Check that the xml version format (schema) is right: IGP, Done
	4. Strand: check the direction of the strand. /IGP, Done/
	5. Build_number: check the number of build provided. If more than one, check that coordinates are diff. /VF,Done/
	6. Build_coor: check that the start and end of the coordinate for each build are different /VF In Progress/
	7. Data_consistency: check that coordinates extracted are numbers /IGP/VF ?/
	8. Exon_coor: check that exons do not overlap /IGP/VF ?/
	9. Exon_number: check that we get the right number of exons /IGP/VF?/
	  
2.Parsing:
Parsing is important or whatever, bla, bla /2-3 lines, VF/

Parsing Features:
/Describe the parsing functions, 1-3 lines per function, IGP + VF as required/

3.Files generation:
Results need to be gathered and collected in a way that allows it further use by other professionals. This could be in 
a human readable way  (e.g. tab separated document, build mapper) or in a script readable way (e.g. bed and cvs file).
... bla, bla 
Information:

Bed file: 
chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display 
of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases
numbered 0-99.

Source: https://genome.ucsc.edu/FAQ/FAQformat#format1

chr22 1000 5000 
chr22 2000 6000


/VF/

File generation features:
/Describe the file generation functions, 1-3 lines per function, IGP + VF as required/

#### STRUCTURE /IGP, Done/

Main program running through different steps: 
1. Initial tests: 
Aim: Run initial tests to check the software and file before execution
Function:initial_tests()

2. Handle xml structure:
Aim: Use 'xml.etree.ElementTreee' to extract information from xml files
Function: (root, up_anno) = handle_xml(data)

3. Getting background information: 
Aim: Get background information about the gene
Function: (schema, lrg_id,  hgnc_id, seq_source, transcript, cs, start_cs, end_cs, strand_cs) = get_background(root)

4. Get information about the different builds: 
Aim: Get build information, including coordinates, chromosome, transcript,and genomic start and end. It will provide "N/A", when protein coordinates are not available
Function: (coord, chro, NC_trans, gstart, gend) = get_build_info(up_anno)

5. Get information about the gene and the strand: 
Aim: Extract gene name (HGVS nomenclature) and tag strand as forward(+) or reverse(-) 
Function:(gene, str_dir) = get_gen_data(data) 

6. Get information about the exon:
Aim: Get information about the exon for the different transcripts, including number of exons, exons coordinates in the LRG system regarding the cdna, transcript and protein.
Function: (list_all_coord, list4bed)= get_exon_data(data, gstart, gend, chro, str_dir)

7. Save exon, transcript and protein coordinates to file
Aim: Creating csv, a tab separate txt file with exon, transcripts and 
protein coordinates and a bed file
Function: output2file(list_all_coord, list4bed)

8. Print disclaimer:
Aim: Print disclaimer statement at the end of program execution
Function: disclaimer()

9. Run of final tests
Aim: Checking for the xml file format/version.
Function: final_tests()

Note: Refer to Developers file for further information


#### SCOPE /VF In Progress/
This software has been tested in Windows and Linux  with a sample of .xml files obtained from the LRG website. This sample contain genes with one or more builds and/or transcripts. 


Genes were selected based on their clinical utility .... the main genes analysed in diagnostic labstories, such as:

/Here we include a list of the genes tested, the gene name and its main impact/

The following genes and their corresponding LRG files were used during testing for reasons outlined below.
LRG_292 = BRCA1: Decisive in the diagnosis of Breast cancer. Used to test script handling of genes that occur on the reverse strand.
LRG_1 = COL1A1:
LRG_9 = SDHD: Used to test script handling of a record that is currently 'Pending Approval' therefore the record is not public yet and subject to change. This was not possible to test as the .xml record for LRGs that are 'Pending Approval' do not contain this phrase within the LRG therefore the only way to know if the LRG is pending approval or has been accepted is to look on the website.
LRG_214 = NF1: Used to test handling of more than one transcript contained within an LRG file.

/ LRG_214 (Gene NF1) Two transcripts LRG_292 (Gene BRCA1) Reverse strand/

#### VERSIONING /VF, IP/

- Feature added: Extracting background information about the gene: v1
- Feature: Making the code modular: v2
- Improvement: Adding examples of initial and final test: v2.1
- Feature: Extracting exon, transcript and protein coordinates: v3
- Feature: Creation of output: csv and tab separate text file: v3.1
- Feature: Dealing with more than one transcript: v.4
- Improvement: Adding disclaimer: v.4.1
- Feature: Creation of a bed file: v.5
- Improvement: Organisation of input and outputs into folders: v.5.1
- Improvement: Splitting of input file name. cleaning of code clean, comments added:v.5.2
- Improvement: Bug resolved for LRG lacking first prot. coordinates in exon 1 
    (e.g. LRG_292 and LRG 62). Comments revised, added and updated v.5.3 
- Feature: Addition of command line arguments. Now ".xml" file can be provided from the command line v.6
- Upgrade: Addition of test for strand direction within function "final_tests" v.6.1
- Feature: v.7 - VF, IP
- Improvement: v.7.1 VF, IP



#### DISCLAIMER /IGP, Done/
Please cite this software as: 'Gomez-Paramio, I. and Fryer, V. (2016), 'lrgext', Software, Faculty of Medicine and Human Sciences, The University of Manchester.' or successor references as defined by the authors.

For further information, refer to COPYRIGHT /IGP, Done/
