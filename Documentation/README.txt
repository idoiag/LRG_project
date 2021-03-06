NOTE: Please refer to https://github.com/idoiag/LRG_project/ to find all the software versions, data sources, outputs and further documentation

#### BACKGROUND ####

Research into variation within the human genome is being undertaken worldwide and is vastly improving our knowledge with regards to benign and disease-causing mutations. However, the rapid pace at which this increase in knowledge is occurring (largely due to the evolution of affordable next generation sequencing) is adding to the already complicated issues of standardising genomic and variation nomenclature. New research leads to better understanding of the size and structure of the human genome, which in turn leads to alterations in numbering of bases, exon/intron boundaries and assessment of variation to find the 'most common' reference genome.

There now exists several genome assembly versions, each subsequent version being slightly different to the previous version. As each version is released, the bioinformatic tools and software used in genetic analysis must be adapted to incorporate the new standard, leading to a time lag between reference genome release and use within the laboratory/clinical setting. As there are several programs used for interpretation of genetic testing, at any one time there may be a combination of more than one reference genome in use within the laboratory/clinic as each development team takes varying lengths of time to incorporate the new standard. It is therefore imperative that, when clinical reports are issued regarding any variants found, the exact location of the the variant(s), the genome assembly (or 'build') and the transcript used to identify the variant are accurately described. This allows any further analysis (or retrospective analysis) to be performed efficiently and accurately as the Scientist has all the relevant information available to them. 
 
Locus Reference Genomic sequences (LRGs) are used as a stable reference for naming variants within a gene and are offered as an alternative (and in the future it is hoped, the standard) method for naming variants. The advantage of the LRG reference system is that an LRG is created for each gene, and once it has been approved by experts for that specific gene in collaboration with both EBI (European Bioinformatics Institute) and NCBI (National Center for Biotechnology Information), the numerical notation will never change for the specific gene. This can be very useful as it means that all notation of a variant will remain the same, allowing easier identification of cases of the same variant within the literature and in the laboratory (particularly important in family investigations that may occur over several generations).

LRG records are available online at www.lrg-sequence.org and can be downloaded and stored locally as an .xml file or FASTA file. LRG records contains a lot of data about each LRG and it is useful for the bioinformatician to be able to access these records and extract a subset of data that they and other clinical scientists within their department may find useful. Developing an xml parser that deals 
specifically with LRG records means that the exact data need can be extracted very quickly and consistently each time the request is made.

Extensible markup language (XML) is a method for transferring data between applications in a structured and standardised format so that the data can be easily converted and read by the local language used (e.g. Python Dictionary or Java HashMap). The data is represented in a hierarchical tree structure and access can be gained to all branches of the tree (i.e. data in those branches) be navigating from the root to the desired branch/information via the 'nodes' in between.

Many laboratories currently store the LRG xml files on a local folder (which will decrease the processing time of the script slightly as it won't have to send requests via the internet) as the LRG xml files are not large. This method is more secure as all the data processing can occur behind the NHS firewall. The disadvantage of this system is that the stored xml files may not be of the most up-to-date information on the website, however this can be minimised by implementing a regular routine for pulling xml files from the site. The script compiles a list of all LRG xml files in the LRGs directory each time it is run before attempting to parse any xml file therefore the script will be using the most up-to-date records that have been saved within the LRGs directory. This means that if old LRG files have been removed, the script will be unable to run them and if new LRG files have been added, the script will automatically recognise them and run them. Therefore, if the LRGs directory is up-to-date, the script will produce up-to-date results.

Parsing xml files using a program that employs defensive programming techniques ensures that the data extracted will be of consistently good quality (i.e. what the Scientist wants) and any errors or unexpected inputs will be dealt with efficiently, with error messages provided to aid the Scientist in using the program effectively. The program script itself will also be easily human readable, with a structure and comments that are clear and concise, which will improve debugging, if it becomes necessary.

#### USAGE/INSTALLATION ####

lrgext does not require installation and can be run from the command line by providing the gene name using the in-built get ops. This script has been customized for most users and does not require knowing in advance the LRG ID, the HGVS gene name (not case sensitive) will be sufficient. Also, it is possible to specify the input folder, containing the LRGs files, the output folder, where outcomes will be saved and the name for the list of genes existing within the input folder. The list will be re-created each time the script is run and thus this script can be used to find out both if an LRG exists within the LRGs directory for a given HGNC gene and its corresponding LRG ID. Flags:
*-g  --gene 
-i --ifile
-o --ofile
-g --lrg

*mandatory. If no arguments are provided, the files will be searched in ./LRGs/, saved in ./Outputs/ and the listing of genes will be �called gene_lrg_lst.csv�.

Example of instruction to the command line: 

$ python lrgext.py -g [HGNC gene name] -i [input folder] -o [output folder] -l [name_desired_for_list_of genes]
$ python lrgext.py -g brAca1
$ python lrgext.py -g BRCA1 -i ./LRGs/ 
$ python lrgext.py -g brca1 -i ./LRGs/ -o ./Outputs/
$ python lrgext.py -g BrcA1 -i ./LRGs/ -o ./Outputs/ -l gene_lrg_lst.csv

N.B. letter case  of the gene name is dealt with within the script therefore any combination of upper/lowercase may be used.

lrgext will extract build, gene, transcript, exon information from a file in LRG xml format producing different output files:
- csv file providing the exon, transcript and protein coordinates separated by comma 
- tab separated txt file providing the exon, transcript and protein coordinates 
- bed file  with the local coordinates and their equivalent genomic coordinates for each exon 
- csv file with information regarding differences between DNA sequences of the LRG and GRCh builds
- tab separated text file with information regarding the GRCh builds

#### FEATURES ####

MAIN FEATURES:

- Versatility: lrgext.py is fully modular making it compatible as a package. Functions can be used as stand-alone by other scripts.
- Efficiency: lrgext provides an efficient extraction of data contained in LRG files using a internationally-recognised standard format("Schema 1.9"). 
- Robustness: due to its modular design, development following "defensive programming" guidance, and extensive testing for bug finding and resolution.
- Compatibility: It has been designed and tested for use in both Linux and Windows environments.
-Flexibility: the incorporation of get opts provide a dynamic way to specify the input file (LRGs directory), the output folder or other features like the name of the LRG files listing file that will be created automatically when the program is run.

1.TESTING: 

Throughout development , testing of the script and the error handling capabilities were frequently and routinely carried out to ensure that the script was extracting the data required, producing files as written and successfully implementing error handling techniques. This testing is very important as there is no way of knowing if the script will work under all the conditions (i.e. different file formats due to different levels of data available for each LRG record) without carrying out extensive testing. 

Testing was performed using python modules and commands specific to these functions (e.g os.path.isfile, os.path.exists, assert) and customized commands using variations of if/elif/else statements and loops where appropriate. Break and exit commands were included where necessary with informative statements being printed to the terminal for the user where results of tests were not sufficient for the programme to continue.

Extensive debugging has been carried out by using several LRG xml files that all have different levels of information associated with them.

== Features:
Testing is divided into 2 sections:

Initial tests:

	1. LRG_ID: check if the given gene has an allocated LRG ID.
	2. LRG file: Check the LRG exists and is a readable file. 
Final tests: 
	3. Schema: Check that the xml version format (schema) is right.
	4. Strand: check the direction of the strand and warn if reverse.
	5. LRG size: check if the LRG size is different compared to both GRCh37 and GRCh38
	6. Gene length: check gene length against reference sequence, provide warning if size difference exists (indicating sequence differences).
	7. Exon_number: check that we get the right number of exons 

  
2.PARSING: 

Being able to parse a file such as an xml or JSON format file is an essential skill for any bioinformatician working in a clinical diagnostic laboratory. This is due to the requirement for laboratories to have accurate and up-to-date information whilst there is an increasing amount of data available. This is especially true in the context of genetic testing due to the impact results can have on the lives of individuals and their families once tested. 

Features:

a) Parsing the name of files contained in directory.
This enables the user to enter the gene name at the command line but the program can continue with the LRG id. Debugging has ensured that if the gene name is not entered accurately or if an xml does not exist in the LRGs directory, the script will exit and an error message informs the user that the gene name entered is not recognised.

b) Parsing xml files. The script uses "ElementTree" (which is available as a Python library) to extract clinically useful information both efficiently and consistently. Since xml files have a defined format (described by its 'schema'), any xml record that uses the same schema will have the same format. Therefore this script can be used for any LRG xml file. 

c) Extract gene name (HGVS nomenclature) and the strand direction. A warning will be provided when genes are on the reverse strand. This serves as a validation method that the data extracted corresponds to the gene initially queried. 
 
d) Extract background information, including the schema, LRG ID, HGVS ID, source of sequence, transcript, coordinate system (the genome assembly) and its start, end and strand coordinates. lrgext.py allows the data extraction of any LRG file (gene) compliant with schema 1.9. (the current and most popular LRG schema), and there will be a warning printed to the terminal if the schema is different.

e) Extract exon, transcript and protein information, including number of exons and their coordinates, transcripts and protein coordinates. These coordinates will be later used to create the comma separated, tab separated and bed files.

e) Extract information about the builds provided within the file, including build number, genomic start and end coordinates, size of the LRG, and size of the gene in different builds. If there is a difference in size of the LRG when compared to any of the reference genomes for which it has information, a warning will be printed to the command line.

f) Extract differences between GRCh builds (which can be any of several reference sequences/builds and more than one comparison may be available) and the LRG itself, with any variations in DNA sequence being highlighted in the diff.csv output file. 

3.REPORTING:

Results need to be collected and stored in a way that allows it to be further used by other professionals. This could be in a human readable way  (e.g. tab separated document or .csv file that can be opened in other software such as MS Excel to improve readability) or in a script readable way (e.g. bed and csv file). The majority of the output from the script will be in the format of .txt, .csv or bed files that can be moved and manipulated as desired, with error handling, warnings and results of checks being printed to the command line for the user to identify if there are any problems immediately.

Features:

a) Creation of a tab separated and a csv file listing coordinates for all exons, transcripts and proteins. 
The file contains 8 columns describing the transcript version, the exon number, exon start and end coordinate, transcript start and end coordinate and protein start and end coordinate. 

Example tab separated file (.txt):
Transcript	Exon	Exon_start	Exon_end	Transcipt_start	Transcript_end	Protein_start	Protein_end
1	1	5001	5229	1	229	1	35

Example comma separated values file (.csv): 
Transcript,Exon,Exon_start,Exon_end,Transcipt_start,Transcript_end,Protein_start,Protein_end
1,1,5001,5229,1,229,1,35

b) Creation of a bed file. Bed files have a specific format, with 3 mandatory columns and 9 optional columns (Source: https://genome.ucsc.edu/FAQ/FAQformat#format1): 
1.chrom - Chromosome number (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
2.chromStart - The starting position of the feature in the chromosome/scaffold. The first base in a chromosome is numbered 0.
3.chromEnd - The end position of the feature in the chromosome/scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.
4.(optional) - strand direction (+/-)
5.(optional) - source, in this case, transcript version

Example:
chrom chromStart chromEnd strand transcript
17 50187097	50211868	-	1

c) Creation of a build comparison tab separated file containing the Build, Chromosome, NC_trans (RefSeq Chromosome accession number),genomic start coordinate, genomic end co-ordinate, LRG start, LRG end, LRG size and gene size.

Example:
Build	chr	NC_trans	geno_start_coord	geno_end_coord	LRG_start	LRG_end	LRG_size	gene_size
GRCh37.p13	17	NC_000017.10	41176312	41370000	1	193689	193688	193688

d) Creation of a sequence differences comma separated values file containing the differences between the DNA sequence of the LRG 
and the DNA sequence of the GRCh builds:
Example:
Build,Diff_type,LRG_start,LRG_end,Geno_coord_start,Geno_coord_end,LRG_allele,Ref_allele
GRCh37.p13,mismatch,859533,859533,32503194,32503194,G,A

#### STRUCTURE ####

The program runs through the following steps: 

1. Create a LRG file list. 
Aim: Create a file listing all LRG files existing on the LRG directory. 
Function: create_repository_file()

2. Initial tests: 
Aim: Tests for the existence of LRG file, their consistency and readability.  
Function: initial_tests()

3. Handle xml structure:
Aim: Use 'xml.etree.ElementTree' to extract information from xml files
Function: (root, up_anno) = handle_xml(data)

4. Get gene and strand information:
Aim: Extract HGVS gene name, and tag strand as forward (+) or reverse (-).
Function: (gene, str_dir) = get_gen_data(data) 

5. Getting background information: 
Aim: Get background information about the gene.
Function: (schema, lrg_id,  hgnc_id, seq_source, transcript, cs, start_cs, end_cs, strand_cs) = get_background(root)

6. Get information about the different builds: 
Aim: Get build information, including coordinates, chromosome, transcript,and genomic start and end. It will provide "N/A", when protein coordinates are not available.
Function: (coord, chro, NC_trans, gstart, gend) = get_build_info(up_anno)

7. Get information about the exon:
Aim: Get information about the exon for the different transcripts, including exons number and coordinates in the LRG system regarding the cdna, transcript and protein.
Function: (list_all_coord, list4bed)= get_exon_data(data, gstart, gend, chro, str_dir)

8. Get differences between LRG sequence and GRCh builds
Aim: Extract information regarding any differences between the LRG sequence and the reference sequence.
Function:(diff_data) = diff_data(data)

9. Save exon, transcript and protein coordinates to file
Aim: Create csv and a tab separated file with exon, transcripts and protein coordinates, plus a bed file.
Function: output2file(list_all_coord, list4bed)

10. Save build information and differences to file
Aim: Create a .csv and a tab separated file highlighting differences between builds including genomic co-ordinates, gene and sequence differences.
Function:diff2file(build_data, diff_data)

11. Run of final tests
Aim: Checking for the xml file format/version.
Function: final_tests()

12. Print disclaimer:
Aim: Print disclaimer statement at the end of program execution.
Function: disclaimer()

13. Main function
Aim: Call all other functions, provide get opts and provide values to variables
Function: main(argv)
Note: Refer to Developers file for further information

#### SCOPE ##### 

This software has been tested in Windows and Linux with a sample dataset of .xml files obtained from the LRG website. This sample dataset contains genes with one or more builds and/or transcripts to test robustness of the script. 

The following genes and associated LRG files were used during testing for the reasons outlined below.

LRG_1 = COL1A1:
Disease-causing variants in this gene are associated with Ehlers-Danlos syndrome, Caffey disease and osteogenesis imperfecta. This LRG record was used to test script handling of a LRG record in which there is a sequence difference between GRCh37 build and the LRG. This means that running this script will produce data within the diff.csv file for GRCh37.

LRG_9 = SDHD: 
This gene is a tumour suppressor and codes for one of four subunits (the others being SDHA, SDHB and SDHC) that comprise the succinate dehydrogenase (SDH) enzyme, important in mitochondrial function. Used to test script handling of a record that is currently 'Pending Approval' therefore the record is not public yet and subject to change. However, during testing it became apparent that this was not possible as the .xml record for LRGs that are 'Pending Approval' do not contain this phrase within the LRG therefore the only way to know if the LRG is pending approval or has been accepted is to look on the website (neither the text 'pending' nor 'approval' is in the xml).

LRG_62 = FOXP3: 
A transcription factor gene found on the X-chromosome, variants in which may cause susceptibility to diabetes mellitus and immunodysregulation (as it is codes a protein that is essential for normal function of T cells). This LRG record was used to test script handling of data from 3 GRCh builds (there is a patch for GRCh37.p13) as well as a difference in size between the LRG and GRCh38 sequence. This will produce a diff.csv with sequence differences from more than one genome assembly, as well as an warning message to the terminal that the LRG size is different to the sequence of GRCh38 (investigation of the build.txt shows they differ by 3 bases between LRG and GRCh38 and 2 bases between LRG and GRCh37).

LRG_85 = MRE11: 
Associated with ataxia-telangiectasia-like disorder. This LRG record was used to general handling of an LRG xml file.

LRG_292 = BRCA1: 
Routinely tested to diagnose breast/ovarian cancer. This LRG record was used to test script handling of an LRG record for a gene that occurs on the reverse strand and in which there are differences between the LRG sequence and the GRCh38 sequence. A warning message is printed to the terminal that states this gene is on the reverse strand and that LRG and GRCh38 differ in size.

LRG_220 = MUTYH:
Associated with colorectal/gastric cancers, particularly and inherited autosomal recessive (AR) form of adenomatous polyposis (benign polyps in the large intestine/rectum that may become cancerous unless treated). This LRG record was used to test script handling of an LRG record for a gene that occurs on the reverse strand, and the LRG size has not changed between builds (the gene size in build 37 is the same as in build 38).

LRG_130 = APC:
Encodes a protein involved in tumour suppression, mutations in  which are linked to several types of cancer (including colorectal cancer, brain tumours and gardner syndrome).This LRG record was used to test script handling of an LRG with more than one transcript (3 in total). This gene has three transcripts, each with a different subset of exons:

t1 = 16 exons (2,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)
t2 = 17 exons (2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)
t3 = 16 exons (1,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)

This record prints 3 lines to the command line, one for each transcript, giving the number of the transcript and the total number of exons in the transcript. Both coord files (.txt and .csv) print the correct exon numbers within each transcript. 

LRG_199 = DMD:
This gene is tested in diagnosing paediatric cardiomyopathies, Duchenne muscular dystrophy and Becker muscular dystrophy. It is one of the largest genes in the human genome and codes for a protein called dystrophin. This LRG record was used to test script handling of variations between LRG sequence and both GRCh37 and GRCh38 sequences.

LRG_214 = NF1: 
Used in diagnosis of neurofibromatosis. This LRG record was used to test script handling of more than one transcript contained within an LRG file.

LRG_517 = RB1:
A tumour suppressor gene that also has a role in apoptosis (programmed cell death) and cell differentiation. Used in the diagnosis of retinoblastoma, also associated with bladder cancer and other form of cancer (lung and breast cancer, osteosarcoma and melanoma). This LRG record was used to test script handling of a gene on the forward strand.

#### FURTHER INFORMATION ####

Information regarding genotype-phenotype relationships and summary of gene functions from OMIM (https://www.omim.org/) and Genetics Home Reference (https://ghr.nlm.nih.gov/).

The LRG records for each gene have been created using a particular DNA sequence, which may or may not be the same as a particular Genome Reference Consortium (GRC) sequence. Once the LRG has been created, the start of the LRG will be at LRG position 1, and the LRG end position will be at the end of the LRG. The size of the LRG will not change once it has been created. However, each time a new GRC human DNA build is released, the LRG entries will be updated to include the genomic start and end coordinates of the gene (which may, and probably will, differ between builds) as well as creating a file that contains any sequence differences between the reference genome and the LRG.

For example, the BRCA1 gene has an LRG (292) that is 193,688 bases long. As it was created using the same sequence as GRCh37, there are no differences between the LRG sequence and the GRCh37 gene sequence. However, the gene coordinates for BRCA1 in GRCh38 have changed such that the gene is now 193,686 bases long. As the LRG sequence does not change, there must be differences between the sequence for the LRG and the gene sequence in GRCh38. These differences are extracted from the LRG xml file into a .csv file called LRG_292_BRCA1_diff.csv. In this file (as in the table on the LRG webpage), the "LRG allele" is the base/sequence present in the LRG and the "Reference allele" is the base/sequence present in the reference sequence being referred to in that instance(in this case, the sequence from GRCh38).

N.B. As more GRCh builds are released, the "reference allele" will change according to which build is being used. The LRG should never change.

#### VERSIONING ####

- Feature added: Extracting background information about the gene v1
- Feature: Making the code modular v2
- Improvement: Adding examples of initial and final test v2.1
- Feature: Extracting exon, transcript and protein coordinates v3
- Feature: Creation of output: csv and tab separate text file v3.1
- Improvement: Bug resolved for LRG lacking first prot. coordinates in exon 1 (e.g. LRG_292 and LRG 62). Comments revised, added and updated v5.3 
- Feature: Addition of command line arguments. Now ".xml" file can be provided from the command line v6
- Upgrade: Addition of test for strand direction within function "final_tests" v6.1
- Feature: Addition of searchable gene as argument v7
- Improvement: Free coded added in functions. Bug linked to Window-Linux environment solved. v7.1 
- Improvement: Build output added v.7.1.1
- Feature: Creation of differences.csv file and exit point if file does not exit v8 
- Improvement: Code simplified, code to print differences build included in new function. v8.1
- Bug fixed: Domain of bed files has been fixed. Usage and comments improved. v8.1.1
- Feature: New tests added v8.2
- Feature: New test added, main function included, structure optimized v8.3
- Feature: Test function added. v8.4
- Improvement: New assert functions added. Output in the command line tidied up v8.4.1
- Improvement: Able to return differences between multiple builds v8.2.1
- Improvement: Differences in multiple builds now captured in diff.csv file (e.g. gene DMD) v8.4.2
- Improvement: Command line tidied up. Hard code eliminated. v8.4.3
- Feature: Get opts added, code cleaned, command line tidied up. Hard code eliminated. v8.5 => lrgext.py

#### DISCLAIMER ####

Please cite this software as: 'Gomez-Paramio, I. and Fryer, V. (2016), 'lrgext', Software, Faculty of Medicine and Human Sciences, The University of Manchester.' or successor references as defined by the authors.

For further information, refer to COPYRIGHT

NOTE: Please refer to https://github.com/idoiag/LRG_project/ to find all the software versions, data sources, outputs and further documentation


