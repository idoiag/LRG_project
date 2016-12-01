/Text within slashes is temp and will be deleted in a final version. IP =In progress. Done = Only minor amendments/addition/updates expected/

##############  README ############

#### BACKGROUND /VF, IP/
/Examples: 
Why LRGs are important, use, challenges,
How can an xml parser help the diagnostic service,
How can bioinformaticians help > great idea to have a parser
Etc… /


#### USAGE/INSTALLATION /IGP, IP, pending on features being added/

lrg does not require installation and can be run from the command line by providing the LRG file as an argument.
lrg will extract build, gene, transcript, exon information from a file in LRG format 
producing different output files:
    - cvs file: provides the exon, transcript and protein coordinates separated by comma /IGP/
    - tab separated txt file: provides the exon, transcrip and protein coordinates separated by tab /IGP/
    - bed file: provides the local coordinates and their equivalent genomic coordinates for each exon /IGP/
    /- build mapper: provides the differences between the builds, highlighting if variants lie in intronic/exonic regions /VF//
	
/get ops code or arguments need still to be added to the program /IGP/. We also need to produce another doc 
highlighting the differences between the builds (VF) /

#### FEATURES /IGP + VF/

/Background:VF 
Example: sell our program highlighting that is robust, efficient and its code could be reused since it was built in modules 
(see features below).
Also, it would be good if we could provide examples of having it built following the principles of defensive programming 
(i.e. using tests, see testing secton) but it will be important to build a test using assert. Ideally we could include a section of how we deal with error handling,


1.Testing: 
Testing is important because bla, bla... / 2-3 lines, VF/

/This is a list of tests we initally suggested. We also proposed to use BRCA as a control. 1-3 lines to describe each
function, IGP  + VF as required/

	1. LRG_ID: check if the given gene has an allocated LRG ID. /VF?/
	2. Check the LRG exists and is a readable file. /IGP, DONE /
	3. Check that the xml version format (schema) is right: IGP, Done

	3. Strand: check the direction of the strand. /IGP, In Progress/
	4. Build_number: check the number of build provided. If more than one, check that coordinates are diff. /VF ?/
	5. Build_coor: check that the start and end of the coordinate for each build are different /VF ?/
	6. Data_consistency: check that coordinates extracted are numbers /IGP/VF ?/
	7. Exon_coor: check that exons do not overlap /IGP/VF ?/
	8. Exon_number: check that we get the right number of exons /IGP/VF?/
	  
2.Parsing:
Parsing is important or whatever, bla, bla /2-3 lines, VF/

Parsing Features:
/Describe the parsing functions, 1-3 lines per function, IGP + VF as required/


3.Files generation:
Results need to be gathered and collected in a way that allows it furhter use by other professionals. This could be in 
a human readable way  (e.g. tab separated document, build mapper) or in a script readable way (e.g. bed and cvs file).
... bla, bla /VF/

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

2. Getting background information: 
Aim: Get background information about the gene
Function: (schema, lrg_id,  hgnc_id, seq_source, transcript, cs, 	start_cs, end_cs, strand_cs) = get_background(root)

3. Get information about the different builds: 
Aim: Get build information, including coordinates, chromosome, transcript,and genomic start and end. It will provide "N/A", when protein coordinates are not available
Function: (coord, chro, NC_trans, gstart, gend) = get_build_info(up_anno)

4. Get information about the gene and the strand: 
Aim: Extract gene name (HGVN nomenclature) and tag strand as forward(+) or reverse(-) 
Function:(gene, str_dir) = get_gen_data(data) 

5. Get information about the exon:
Aim: Get information about the exon for the different transcripts, including number of exons, exons coordinates in the LRG system regarding the cdna, transcript and protein.
Function: (list_all_coord, list4bed)= get_exon_data(data, gstart, gend, chro, str_dir)

6. Save exon, transcript and protein coordinates to file
Aim: Creating csv, a tab separate txt file with exon, transcripts and 
protein coordinates and a bed file
Function: output2file(list_all_coord, list4bed)

7. Print disclaimer:
Aim: Print disclaimer statement at the end of program execution
Function: disclaimer()

8. Run of final tests
Aim: Checking for the xml file format/version.
Function: final_tests()

Note: Refer to Developers file for further information

#### SCOPE /VF/
This software has been tested with a sample of .xml files obtained from the LRG website. 
We have tried to include the main genes analysed in diagnostic labstories, such as:

/Here we include a list of the genes tested, the gene name and its main impact/

LRG_292 – BRCA1 – Decisive in the diagnosis of Breast cancer
LRG_ ….


/N.B. LRG.xml files included for testing purposes include: LRG_1 (Gene COL1A1) LRG_9 (Gene SDHD) Pending approval, not public yet LRG_214 (Gene NF1) Two transcripts LRG_292 (Gene BRCA1) Reverse strand/


#### VERSIONING /IGP, Done/

- Extracting background information about the gene: v1
- Making the code modular: v2
- Adding examples of initial and final test: v2.1
- Extracting exon, transcript and protein coordinates: v3
- Output: it creates a csv and tab separate text file: v3.1
- Dealing with more than one transcript: v.4
- Adding disclaimer: v.4.1
- Creation of a bed file: v.5
- Organisation of input and outputs into folders: v.5.1
- Automatic creation of output file name. Code clean, comments added:v.5.2
- Bug resolved for LRG lacking first prot. coordinates in exon 1 
    (e.g. LRG_292 and LRG 62). Comments revised, added and updated v.5.3 


#### DISCLAIMER /IGP, Done/
Please cite this software as: 'Gomez-Paramio, I. and Fryer, V. (2016), 'lrgext', Software, Faculty of Medicine and Human Sciences, The University of Manchester.' or successor references as defined by the authors.

For further information, refer to COPYRIGHT /IGP, Done/


     
 

