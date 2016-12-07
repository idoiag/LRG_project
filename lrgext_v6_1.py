# -*- coding: utf-8 -*-
# Authors: Idoia Gomez-Paramio and Verity Fryer 2016
# Usage: python lrgext_v6_1.py LRG_ID (e.g. LRG_292)

import xml.etree.ElementTree as ET
import os.path
import sys


"""
lrgext extracts build, gene, transcript, exon information from a file in LRG format
producing different output files:
    - .csv file
    - tab separated .txt file
    - bed file
"""

"""
Capturing file and initializing variables
"""

script = sys.argv[0]
LRG = sys.argv[1]
#LRG = 'LRG_62' for example
path = './LRGs/'
data = path + LRG + '.xml'
# Add error checking to detect if LRG doesn't exist
#    print("No .xml file exists for the name specified")

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
    print('Gene: ', gene)

    # Determining the strand direction(e.g. 517 is "+")
    for annotation in root.findall('./updatable_annotation/annotation_set[@type="lrg"]/mapping[@type="main_assembly"]/mapping_span'):
        str_dir = annotation.attrib['strand']

        # marking forward strand as "+" and negative as "-"
        if (str_dir == "1"):
            str_dir = "+"
        else:
            str_dir = "-"
        print (str_dir)

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

        print ("Schema version: " + schema)
        print ("LRG  ID: " + lrg_id + "\n", hgnc_id, seq_source, transcript, cs, start_cs, end_cs, strand_cs)
        return (schema, lrg_id, hgnc_id, seq_source, transcript, cs, start_cs, end_cs, strand_cs)

"""
Get build information, including coordinates, chromosome, transcript, genomic start and genomic end.
It will provide "N/A", when protein coordinates are not available
"""
def get_build_info(up_anno):

    for annotation in up_anno[1].findall('mapping'):
        build = annotation.get('coord_system')
        chro = annotation.get('other_name')
        NC_trans = annotation.get('other_id')
        gstart = annotation.get('other_start')
        gend = annotation.get('other_end')
        print(build, chro, NC_trans, gstart, gend)
    return(build, chro, NC_trans, gstart, gend)


####VF CODE###
        # Not sure if the differences between builds could be included here or in a different function. 
        # Also to capture if difference in sequences occur in the intron/exon
###

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
        print( "Transcript number: ", trans_number, 'Exon count: ', len(exon_lst))

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
                    print ("Problem extracting exon information")

            # Create list of coordinates
            list4bed.append([chro, g_start_ex, g_end_ex, str_dir, str(trans_number)])
            list_all_coord.append([str(trans_number), ex_num, start_ex, end_ex, start_ex_tr, end_ex_tr, start_ex_pt,end_ex_pt]  )

    # Preparing lists to be print in columns
    for group in list_all_coord:
        print ("\t".join(group) + "\n")
    return (list_all_coord, list4bed)

"""
Creating .csv, a comma separated text file with exon, transcripts, protein coordinates and a bed file"""

def output2file(list_all_coord, list4bed):
    # Open file in read/write format. If one doesn't exist, it  will create a new one
    db = open("./Outputs/LRG_coord.txt","w")
    db_csv = open("./Outputs/LRG_coord.csv","w")
    bed = open ("./Outputs/LRG_bed", "w")

    # Add headers to files
    headings = ["transcript","exon", "ex_start", "ex_end", "tr_start", "tr_end", "pt_start", "pt_end"]
    bed_headings = ["chr", "start", "end", "strand", "transcript" ]

    # Writing tab separated text file
    db.write("\t".join(headings) + "\n") # writing headings
    for group in list_all_coord:
        db.write ("\t".join(group) + "\n") # writing coordinates

    # Writing csv file
    db_csv.write(",".join(headings) + "\n") # writing headings
    for group in list_all_coord:
        db_csv.write (",".join(group) + "\n") # writing coordinates

    # Writing bed file
    bed.write("\t".join(bed_headings) + "\n") # writing headings
    for group in list4bed:
        bed.write ("\t".join(group) + "\n")

    # Closing files
    db.close()
    db_csv.close()
    bed.close()
    return

def disclaimer():
# This disclaimer is in triple quotes therefore won't print - is this because i needs to be in double quotes or because it needs to be amended first?
    print ("""Please cite this software as: 'Gomez-Paramio, I. and Fryer, V. (2016), 'lrgext', Software, 
    Faculty of Medicine and Human Sciences, The University of Manchester.' or successor 
    references as defined by the authors.\n""")
    return

##### TESTING #####

def initial_tests():
    """ Run initial tests to check the software and file before execution""" 
    if (os.path.isfile(data) == False) :
        print ("Data is not a readable file")

    if (os.path.exists(data) == False):
        print ("Data does not exist")

#### VF's code ####
    # check file is in xml format. If not, return error message "Not an xml file"
    # add try, except to close program if no LRG exists
    return

"""
Tests run at the end of the program
"""
def final_tests():

    """ Checking for xml format """
    if schema != "1.9":
        print ("""lrgext supports LRG xmls built using schema 1.9. Please be aware of data incongruencies""")

    # Check the strand direction and warn if in reverse.
    # Use LRG_571 for a forward strand example
    if str_dir == "-":
        print ("Note: reverse strand!")

    return


def builds():

    for annotation in up_anno[1].findall('mapping'):
        build = annotation.get('coord_system')

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

            print(build, NC_trans, gstart, gend, lrg_start, lrg_end, lrg_size_37)

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

            print(build, NC_trans, gstart, gend, lrg_start, lrg_end, lrg_size_38)

        # if any build other than 37 or 38 is present, this script will need to be modified
        else:
            print("This is not an expected Human Reference Genome Build. Please check .xml file")
            break
    if lrg_size_38 == lrg_size_37:
        print("LRG size is the same between GRCh37 and GRCh38")
    else:
        print("WARNING: LRG sizes differ between GRCh37 and GRCh38")


#### VF's code
    # if transcript number differ between builds 37 and 38:
    # print("Do transcripts match between builds? ", True/False)
    return

### End of testing ###

#### MAIN information ####

initial_tests()   # 1
(root, up_anno) = handle_xml(data)  # 2
(schema, lrg_id,  hgnc_id, seq_source, transcript, cs, start_cs, end_cs, strand_cs) = get_background(root) # 3
(build, chro, NC_trans, gstart, gend) = get_build_info(up_anno) # 4
(gene, str_dir) = get_gen_data(data) # 5
(list_all_coord, list4bed) = get_exon_data(data, gstart, gend, chro, str_dir) # 6
output2file(list_all_coord, list4bed) # 7
disclaimer() # 8
final_tests() # 9
builds()

### End of Main ###
