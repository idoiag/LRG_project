# -*- coding: utf-8 -*-
# Authors: Idoia Gomez-Paramio and Verity Fryer 2016
# Usage: python lrgext_v7.1.py enter_gene (e.g. BRCA1, maybe be uppercase or lowercase or a combination of both)

"""
lrgext extracts build, gene, transcript, exon information from a file in LRG format
producing different output files:
    - .csv file
    - tab separated .txt file
    - bed file
"""

import xml.etree.ElementTree as ET, os.path, sys, csv


"""
Capturing file and initializing variables
"""

script = sys.argv[0]
#LRG = sys.argv[1]
path = './LRGs/'

enter_gene = sys.argv[1].upper()

# parse all .xml files in LRGs directory and create a list in .csv format of all LRG file names within the directory and their associated gene names
with open('gene_lrg_lst.csv','w') as f:
    csvfileWriter = csv.writer(f)
    for filename in os.listdir(path):
        if not filename.endswith('.xml'): continue
        fullname = os.path.join(path, filename)
        tree = ET.parse(fullname)
        root = tree.getroot()
        gene = root.find('updatable_annotation/annotation_set/lrg_locus').text

        csvfileWriter.writerow([gene, filename])
    f.close()

##### TESTING #####

def initial_tests():
    """ Run initial tests to check the software and file before execution""" 
    if (os.path.isfile(data) == False) :
        print ("\nData is not a readable file")

    if (os.path.exists(data) == False):
        print ("\nData does not exist")

    return

# read the .csv file of LRGs to check an LRG exists for the gene name entered
with open('gene_lrg_lst.csv','r') as f:
    csvfileReader = csv.reader(f)
    # Iterate through .csv file and check first entry in each row for gene name
    for row in csvfileReader:
        found = 'No'
        if enter_gene in row[0]:
            found = 'Yes'
            LRG_xml = row[1]
            # if the gene name entered is in the first entry for a row, exit the loop
    if found == 'No':
        print("\nNo LRG.xml can be found for gene entered.\nCheck that gene has an associated LRG at www.lrg-sequence.org")
        break
data = path + LRG_xml

#A further improvement would be the addtion of get -ops
#https://www.tutorialspoint.com/python/python_command_line_arguments.htm
#getopt.getopt(LRG, -l:, [long_options])

# Parsing .xml file with 'xml.etree.ElementTree' 

def handle_xml(data):
    tree = ET.parse(data)
    root = tree.getroot()
    up_anno = tree.getroot()[1]
    return(root, up_anno)

"""
# Extract gene name (HGVS nomenclature) and tag strand as forward(+) or reverse(-)
"""
def get_gen_data(data):

    # Extracting name of gene
    gene = root.find('updatable_annotation/annotation_set/lrg_locus').text
    print('\nGene: ', gene)

    # Determining the strand direction(e.g. 517 is "+")
    for annotation in root.findall('./updatable_annotation/annotation_set[@type="lrg"]/mapping[@type="main_assembly"]/mapping_span'):
        str_dir = annotation.attrib['strand']

        # marking forward strand as "+" and negative as "-"
        if (str_dir == "1"):
            str_dir = "+"
        else:
            str_dir = "-"

    return (gene, str_dir)

"""
Get gene background information 
"""
def get_background(root):

    for lrg in root.findall ("."):
        schema = lrg.get('schema_version')

    for fixed in root.findall("./fixed_annotation"):
        lrg_id = fixed.find('id').text
        hgnc_id = fixed.find ('hgnc_id').text
        seq_source = fixed.find ('sequence_source').text

        for transcript in root.findall("./fixed_annotation/transcript"):
            transcript = transcript.get('name')

            path_fix_coor = "./fixed_annotation/transcript/coordinates"
            for coordinates in root.findall(path_fix_coor):
                cs = coordinates.get('coord_system')
                start_cs = coordinates.get('start')
                end_cs = coordinates.get('end')
                strand_cs = coordinates.get('strand')

        print ("\nSchema version: " + schema)
        print ("\nLRG ID: " + lrg_id)
        return (schema, lrg_id, hgnc_id, seq_source, transcript, cs, start_cs, end_cs, strand_cs)

# Get build information, including coordinates, chromosome, transcript, genomic start and genomic end.
# It will provide "N/A", when protein coordinates are not available

def get_build_info(up_anno):

    for annotation in up_anno[1].findall('mapping'):
        build = annotation.get('coord_system')
        chro = annotation.get('other_name')

        # collect info on build 37
        if build.startswith('GRCh37'):
            NC_trans = annotation.get('other_id')
            gstart = annotation.get('other_start')
            gend = annotation.get('other_end')

            # determine start and end of LRG
            for lrg in up_anno[1].findall('mapping/mapping_span'):
                lrg_start = int(lrg.get('lrg_start'))
                lrg_end = int(lrg.get('lrg_end'))
                if lrg_end > lrg_start:
                   lrg_size_37 = (1+(lrg_end - lrg_start))
                else:
                   lrg_size_37 = (1+(lrg_start - lrg_end))

            print('\n' + build, NC_trans, gstart, gend, lrg_start, lrg_end)

        # otherwise collect info on build 38
        elif build.startswith('GRCh38'):
            NC_trans = annotation.get('other_id')
            gstart = annotation.get('other_start')
            gend = annotation.get('other_end')

            # determine start and end of LRG
            for lrg in up_anno[1].findall('mapping/mapping_span'):
                lrg_start = int(lrg.get('lrg_start'))
                lrg_end = int(lrg.get('lrg_end'))
                if lrg_end > lrg_start:
                   lrg_size_38 = (1+(lrg_end - lrg_start))
                else:
                   lrg_size_38 = (1+(lrg_start - lrg_end))

            print('\n' + build, NC_trans, gstart, gend, lrg_start, lrg_end)

        # if any build other than 37 or 38 is present, this script will need to be modified
        else:
            print("\nThis is not an expected Human Reference Genome Build. Please check .xml file")
            break
    if lrg_size_38 == lrg_size_37:
        print("\nLRG size is the same between GRCh37 and GRCh38")
    else:
        print("WARNING: LRG sizes differ between GRCh37 and GRCh38")
    return(build, chro, NC_trans, gstart, gend, lrg_start, lrg_end)

#def is_int(gstart, gend):
#    assert str(gstart).isdigit()
#    assert str(gend).isdigit()
#    return(gstart_int,gend_int)

"""
Get information about the exon for the different transcripts, including number of exons, exons coordinates in the LRG system regarding the cdna, transcript and protein
"""
def get_exon_data(data, gstart, gend, chro, str_dir):

    trans_number = 0
    list_all_coord, list4bed = [], []

    # Loop to extract information when there is more than 1 transcript (e.g. LRG_214 has two transcripts)
    for transcripts in root.findall('./fixed_annotation/transcript'):
        trans_number += 1

        exon_lst = root.findall('./fixed_annotation/transcript/exon')
        print("\nNumber of transcripts: ", trans_number)
        print('\nExon count: ', len(exon_lst))

        for exons in transcripts.findall('exon'):
            ex_num = exons.get('label')

            # Extracting transcript, protein and exon coordinates (in this order)
            # Print "N/A" if protein coordinates are not available (e.g. LRG_292)
            start_ex_pt = "N/A"
            end_ex_pt = "N/A"

            for coord in exons.findall('coordinates'):
                if (coord.get('coord_system').find("t")!=-1):
                    start_ex_tr = coord.get('start')
                    end_ex_tr = coord.get('end')

                elif (coord.get('coord_system').find("p")!=-1):
                    start_ex_pt = coord.get('start')
                    end_ex_pt = coord.get('end')

                elif (coord.get('coord_system').find("p")==-1) or (coord.get('coord_system').find("t" )==-1):
                    start_ex = coord.get('start')
                    end_ex = coord.get('end')
                    g_start_ex = start_ex + gstart
                    g_end_ex = end_ex + gend

                else:
                    print ("\nProblem extracting exon information")

            # Create list of coordinates
            list4bed.append([chro, g_start_ex, g_end_ex, str_dir, str(trans_number)])
            list_all_coord.append([str(trans_number), ex_num, start_ex, end_ex, start_ex_tr, end_ex_tr, start_ex_pt,end_ex_pt]  )

    # Preparing lists to be printed in columns
    for group in list_all_coord:
        pass
        # print ("\t".join(group) + "\n")
    return (list_all_coord, list4bed)

# Parse xml and retrieve data for all sequence differences between build 37 and 38

def diff_data(data):
    for diff in root.findall('./updatable_annotation/annotation_set[@type="lrg"]/mapping[@type="main_assembly"]/mapping_span/diff'):
        diff_type = diff.get('type')
        diff_lrg_start = diff.get('lrg_start')
        diff_lrg_end = diff.get('lrg_end')
        diff_gen_start = diff.get('other_start')
        diff_gen_end = diff.get('other_end')
        lrg_seq_base = diff.get('lrg_sequence')
        other_seq_base = diff.get('other_sequence')
        print(diff_type,diff_lrg_start,diff_lrg_end,diff_gen_start,diff_gen_end,lrg_seq_base,other_seq_base)
    return(diff_type,diff_lrg_start,diff_lrg_end,diff_gen_start,diff_gen_end,lrg_seq_base,other_seq_base)

"""
Creating .csv, a comma separated text file with exon, transcripts, protein coordinates and a bed file"""

def output2file(list_all_coord, list4bed):
    # Open file in read/write format. If one doesn't exist, it  will create a new one
    db = open("./Outputs/" + lrg_id + "_" + enter_gene + "_coord" + ".txt","w")
    db_csv = open("./Outputs/" + lrg_id + "_" + enter_gene + "_coord" + ".csv","w")
    bed = open ("./Outputs/" + lrg_id + "_" + enter_gene + "_bed", "w")
#    db_build = open("./Outputs/" + lrg_id + "_" + enter_gene + "_build" + ".txt","w")

    # Add headers to files
    headings = ["Transcript","Exon", "Exon_start", "Exon_end", "Transcipt_start", "Transcript_end", "Protein_start", "Protein_end"]
    bed_headings = ["Chromosome", "Start", "End", "Strand", "Transcript" ]
#    build_headings = ["Build", "LRG_start", "LRG_end"]

    # Writing tab separated text file
    db.write("\t".join(headings) + "\n") # writing headings
    for group in list_all_coord:
        db.write ("\t".join(group) + "\n") # writing coordinates

    # Writing comma-seperated values (.csv) file
    db_csv.write(",".join(headings) + "\n") # writing headings
    for group in list_all_coord:
        db_csv.write (",".join(group) + "\n") # writing coordinates

    # Writing tab-delimited .bed file
    bed.write("\t".join(bed_headings) + "\n") # writing headings
    for group in list4bed:
        bed.write ("\t".join(group) + "\n")

    # Writing build tab-delimited text file
#    db_build.write("\t".join(build_headings) + "\n")
#    for group in def_builds:
#        db_build.write("\t".join(grou) + "\n")

    # Closing files
    db.close()
    db_csv.close()
    bed.close()
#    db_build.close()
    # Print to screen creation of output files has been successful
    print ("\nOutput files have been saved in /Outputs!")

    return

def disclaimer():
    print ("""\nPlease cite this software as: Gomez-Paramio, I. and Fryer, V. (2016), 'lrgext', Software, 
    Faculty of Medicine and Human Sciences, The University of Manchester.' or successor 
    references as defined by the authors.\n""")
    return


"""
Tests run at the end of the program
"""
def final_tests():

    """ Checking for xml format """
    if schema != "1.9":
        print ('\nlrgext supports LRG xmls built using schema 1.9. Please be aware of data incongruencies')

    # Check the strand direction and warn if in reverse.
    # Use LRG_571 for a forward strand example
    if str_dir == "-":
        print("\nN.B. This LRG is on the REVERSE strand")
    elif str_dir == '+':
        print("\nN.B. This LRG is on the FORWARD strand")
    return

### End of testing ###

#### MAIN information ####

initial_tests()   # 1
#(gstart_int,gend_int) = is_int(gstart,gend)
(root, up_anno) = handle_xml(data)  # 2
(schema, lrg_id,  hgnc_id, seq_source, transcript, cs, start_cs, end_cs, strand_cs) = get_background(root) # 3
(gene, str_dir) = get_gen_data(data) # 4
(build, chro, NC_trans, gstart, gend, lrg_start, lrg_end) = get_build_info(up_anno) # 5
(list_all_coord, list4bed) = get_exon_data(data, gstart, gend, chro, str_dir) # 6
(diff_type,diff_lrg_start,diff_lrg_end,diff_gen_start,diff_gen_end,lrg_seq_base,other_seq_base) = diff_data(data)
output2file(list_all_coord, list4bed) # 7
final_tests() # 8
disclaimer() # 9

### End of Main ###
