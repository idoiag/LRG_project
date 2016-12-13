/Text within slashes is temp and will be deleted in a final version. IP = In progress. Done = Only minor amendments/addition/updates expected/

##############  README ############

#### BACKGROUND 
Research into variation within the human genome is being undertaken worldwide and is vastly improving our 
knowledge with regards to benign and disease-causing mutations at a rapid pace. However, this exponential 
increase in knowledge (in part due to the evolution of next generation sequencing) is adding to the 
already complicated issues of standardising genomic co-ordinates and variation nomenclature. New research 
leads to better understanding of the genomic map, which in turn leads to alterations in numbering of bases, 
exon/intron boundaries and assessment of variation to find the 'most common' reference genome. There now 
exists several genome assembly versions, each of which has slightly different numbering depending on the 
sequence being investigated. Therefore, when clinical reports are issued regarding any variants found, it 
is very important for the exact location of the the variant(s) to be clearly identified, the genome 
assembly, or 'build' and the transcript used identify the variant. This allows any further analysis 
(or indeed retrospective analysis) to be performed efficiently and accurately as the Scientist has all 
the relevant information to them. However, changes in builds and different transcipts used still causes 
confusion and it can take some time to unravel the nomenclature used between builds.
 
Locus Reference Genomic sequences (LRGs) are used as a stable reference for naming variants within a gene 
and are offered as an alternative (and in the future it is hoped, the standard) method for naming variants. 
The advantage of the LRG reference system is that an LRG is created for each gene, and once it has been 
approved by experts for that specific gene in collaboration with both EBI (European Bioinformatics 
Institute) and NCBI (National Center for Biotechnology Information), the numbering system will never 
change, eliminating the need for clinical reports to contain so much confusing information regarding 
builds and nomenclature and improving consistency and quality of variant reporting. This will also mean 
that each report generated for a specific gene will always have the same LRG associated with it, which can 
be useful in many situations, including family investigations that may occur over several generations.

LRG records are available online at www.lrg-sequence-org and can be downloaded and stored locally as an 
.xml file or FASTA file. LRG records contains a lot of data about each LRG and it is useful for the 
bioinformatician to be able to access these records and extract a subset of data that they and other 
clinical scientists within their department may find useful. Developing an xml parser that deals 
specifically with LRG records means that the exact data need can be extracted very quickly and 
consistently each time the request is made.

Extensible markup anguage (XML) is a method for transferring data between applications in a structured and 
standised format so that the data can be easily converted and read by the local language used (e.g. Python 
Dictionary of Java HashMap). The data is represented in a hierarchical tree structure and access can be 
gained to all branches of the tree (i.e. data in those branches) be navigating from the root to the 
desired branch/information via the 'nodes' in between.

Many laboratories currently store the LRG xml files on a local folder (which will decrease the processing 
time of the script slightly as it won't have to send requests via the internet). This method is also more 
secure as all the data processing can occur behind the NHS firewall. The disadvantage of this system is 
that the stored .xml files may not be of the most up-to-date information on the website, but this can be 
minimised by implemeting a regular routine for pulling .xml files from the site. The script itself 
compiles a list of all LRG xml files in the LRGs directory at the start of the script each time it is run 
therefore the script will be using the most up-to-date records that have been saved within the LRGs 
directory. This means that if additional or amended xml files are added to the directory, the script will 
use these and it won't fall over because it doesn't recognise a new LRG name or //?// if one is missing 
that had previously been processed.

#### USAGE/INSTALLATION /IGP, IP, pending on features being added/

lrgext does not require installation and can be run from the command line by providing the gene name (not case 
sensitive)as an argument (default extension is .xml). This script was originally using the LRG name as the argument, 
but as most users of the script will be most likely to know the HGNC gene name before the LRG ID, this 
script can be used to find out a) if an LRG exists within the LRGs directory for a given HGNC gene and b) 
what the LRG ID is for a given gene name. Example: 

$ python lrgext [HGNC gene name]

lrgext will extract build, gene, transcript, exon information from a file in LRG format 
producing different output files:
    - cvs file: provides the exon, transcript and protein coordinates separated by comma 
    - tab separated txt file: provides the exon, transcrip and protein coordinates separated by tab 
    - bed file: provides the local coordinates and their equivalent genomic coordinates for each exon 
    /- build mapper: provides the differences between the builds, highlighting if variants lie in 
intronic/exonic regions /VF//
	
/We  need to produce another doc highlighting the differences between the builds (VF) /

#### FEATURES /IGP + VF/

/Background:VF 
Example: sell our program highlighting that it is robust, efficient and its code could be reused since 
it was built in modules (see features below).
Also, it would be good if we could provide examples of having it built following the principles of 
defensive programming (i.e. using tests, see testing secton) but it will be important to build a test 
using assert. Ideally we could include a section of how we deal with error handling,

1.TESTING: 
Testing is important because bla, bla... / 2-3 lines, VF/

/This is a list of tests we initally suggested. We also proposed to use BRCA as a control. 1-3 lines to 
describe each function, IGP  + VF as required/

== Features:

	1. LRG_ID: check if the given gene has an allocated LRG ID. 
	2. LRG file: Check the LRG exists and is a readable file. 
	3. Schema: Check that the xml version format (schema) is right.
	4. Strand: check the direction of the strand and warn if reverse.
	5. Build_number: check the number of builds. If more than one, check that coordinates are diff. 
	6. Build_coor: check that the start and end of the coordinate for each build are different /VF In Progress/
	7. Data_consistency: check that coordinates extracted are numbers /IGP/VF ?/
	8. Exon_coor: check that exons do not overlap /IGP/VF ?/
	9. Exon_number: check that we get the right number of exons /IGP/VF?/
	  
2.PARSING: 
Being able to parse an xml file is a today's need for any bioinformatician working in a clinical 
diagnostic laboratory. The reason behind is increasing dependency of lab on accurate and updated information
and the increasing amount of data available over the time. 
This is especially true in the context of genetic testing due to the impact results can have on 
the lives of individuals and their families once tested. 

== Features:
a) Parsing xml files.
The script uses "ElementTree (which is available as a Python library) to extract clinically useful 
information both efficently and consistently. Since xml files have a defined format described by its 'schema'), 
any xml record that uses the same schema will have the same format. 

b) Extract gene name (HGVS nomenclature) and the strand direction. A warning will be provided when genes are 
on the reverse strand. This serve as validation that the data extracted corresponds to the gene initally
queried. 
 
c) Extract background information, including including the schema, LRG ID, HGVS ID, source of sequence, 
transcript, the coordinate system and its start, end and strand coordinates. 
lrgext allows the data extraction of any LRG file (gene) compliant schema 1.9. (current and most popular 
LRG schema), and will worn if the schema is different.

///Are we going to make any use of the coordinate system ??? //

d) Extract exon, transcript and protein information, including number of exons and their coordinates, 
transcripts and protein coordinates. These coordinates will be later used to create the cvs, tab separated and 
bed files.

e) Extract information about the builds provided within the file, including... /VF /

f) Extract differences between builds /VF /

3.REPORTING:
Results need to be gathered and collected in a way that allows it further use by other professionals. 
This could be in a human readable way  (e.g. tab separated document, build mapper) or in a script readable 
way (e.g. bed and cvs file). ... bla, bla  /VF, 3-4 lines will be enough/

== Features: /IGP + VF as required/

a) Creation of a tab separated and a cvs file listing coordinates for all exons, transcripts and proteins. 
The file contains 9 columns describing the transcrip version, the exon number, exon start and end coordinate,
transcript start and end coordinate and protein start and end coordinate. Example:

- Tab separated file (.txt):
Transcript	Exon	Exon_start	Exon_end	Transcipt_start	Transcript_end	Protein_start	Protein_end
1	1	5001	5229	1	229	1	35

- cvs file:
Transcript,Exon,Exon_start,Exon_end,Transcipt_start,Transcript_end,Protein_start,Protein_end
1,1,5001,5229,1,229,1,35

b) Creation of a bed file. Bed files have a specific format, with 3 mandatory columns and 9 optional columns: 
Source: https://genome.ucsc.edu/FAQ/FAQformat#format1 /IGP, Need to make this bit shorter/

1.chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
2.chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a 
chromosome is numbered 0.
3.chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not 
included in the display of the feature. For example, the first 100 bases of a chromosome are defined as 
chromStart=0, chromEnd=100, and span the bases numbered 0-99.
4.(optinal) - strand direction
5.(optional) - source, in this case, transcript version

Example:
chrom chromStart chromEnd strand transcript
17 50187097	50211868	-	1

c) Creation of a build mapper /or whatever techny name you think about :) VF /

d) ... /VF/

#### STRUCTURE /IGP, Done/

Main program running through different steps: 

1. Initial tests: 
Aim: Run initial tests to check for the existence of an LRG file for the given gene name.
It also checks the consistency and readabilyt of the file. 
Function:initial_tests()

2. Create a LRG file list. 
Aim: Create a file listing all LRG files existing on the LRG directory. 
Function:create_repository_file()

3. Handle xml structure:
Aim: Use 'xml.etree.ElementTreee' to extract information from xml files
Function: (root, up_anno) = handle_xml(data)

4. Getting background information: 
Aim: Get background information about the gene
Function: (schema, lrg_id,  hgnc_id, seq_source, transcript, cs, start_cs, end_cs, strand_cs) = 
get_background(root)

5. Get gene and strand information:
Aim:Extract HGVS gene name, and tag strand as forward (+) or reverse (-) 
Function: (gene, str_dir) = get_gen_data(data) 

6. Get information about the different builds: 
Aim: Get build information, including coordinates, chromosome, transcript,and genomic start and end. It will provide "N/A", when protein coordinates are not available
Function: (coord, chro, NC_trans, gstart, gend) = get_build_info(up_anno)

7. Get differences between the builds
Aim: /VF, as short as possible, no more than 1 line =80 characters)
Function:(diff_data) = diff_data(data)

8. Get information about the exon:
Aim: Get information about the exon for the different transcripts, including number of exons, exons coordinates in the LRG system regarding the cdna, transcript and protein.
Function: (list_all_coord, list4bed)= get_exon_data(data, gstart, gend, chro, str_dir)

7. Save exon, transcript and protein coordinates to file
Aim: Creating csv, a tab separate txt file with exon, transcripts and protein coordinates and a bed file
Function: output2file(list_all_coord, list4bed)

8. Save build differences to file
Aim: Create a comma and a tab separated file highlighting difference between builds
Function:diff2file(build_data, diff_data)

9. Print disclaimer:
Aim: Print disclaimer statement at the end of program execution
Function: disclaimer()

10. Run of final tests
Aim: Checking for the xml file format/version.
Function: final_tests()

Note: Refer to Developers file for further information


#### SCOPE /VF In Progress/
This software has been tested in Windows and Linux  with a sample of .xml files obtained from the LRG 
website. This sample contain genes with one or more builds and/or transcripts. 

Genes were selected based on their clinical utility .... the main genes analysed in diagnostic labstories, 
such as:

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
- Feature: /VF, IP/ v.7 
- Improvement:  /VF, IP/v.7.1 
- Improvement: Free coded added in functions. Bug linked to Window-Linux environment solved. v.7.2
- Feature: /VF/ v.8
- Improvement: Code simplified, code to print differences build included in new function. v.8.1




#### DISCLAIMER /IGP, Done/
Please cite this software as: 'Gomez-Paramio, I. and Fryer, V. (2016), 'lrgext', Software, Faculty of Medicine and Human Sciences, The University of Manchester.' or successor references as defined by the authors.

For further information, refer to COPYRIGHT /IGP, Done/
