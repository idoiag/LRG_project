
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

#script = sys.argv[0]
#LRG = sys.argv[1]
path = './LRGs/'
opath = './Outputs/'

#enter_gene = sys.argv[1].upper()
enter_gene = "APC" # For testing purposes

def create_repository_file():
    """ Parse all .xml files in LRGs directory ("LRGs") and  creates a csv file 
    ("gene_lrg_lst.csv") listing all existing .xml files within the directory. 
    The file contains a column with the name of the gene and another
    with its associated .xml file (named with its LRG identifier)
    """
    with open('gene_lrg_lst.csv','w',  newline='') as f:
        csvfileWriter = csv.writer(f)
        for filename in os.listdir(path):
            if not filename.endswith('.xml'): continue
            fullname = os.path.join(path, filename)
            tree = ET.parse(fullname)
            root = tree.getroot()
            gene = root.find('updatable_annotation/annotation_set/lrg_locus').text

            csvfileWriter.writerow([gene, filename])

        f.close()
    return

##### TESTING #####

def initial_tests():
    """ Run initial tests to check the software and file before execution: 
    1st. The LRG file for the povided gene will be searched in a list containing 
    all the LRG files downloaded. If the file exist, the LRG ID will be captured.
    If not, a warning will be displayed. 
    2nd. data consistency within the file will be checked.
    """ 

    # Checking if file  exist and capture name of file 
    with open('gene_lrg_lst.csv','r') as f:
        csvfileReader = csv.reader(f)
        # Iterate through .csv file and check first entry in each row for gene name
        for row in csvfileReader:
            found = 'No'
            if enter_gene in row[0]:
                found = 'Yes'
                LRG_xml = row[1]
                # if the gene name entered is in the first entry for a row, exit the loop
                break

        if found == 'No':
            # If no LRG.xml found for the gene enetered, display an error message and exit the script immediately
            print("\nNo LRG.xml can be found for gene entered.\nCheck that gene has an associated LRG at www.lrg-sequence.org\n")
            exit()

    data = path + LRG_xml

    # Check data consistency within file 
    if (os.path.isfile(data) == False) :
        print ("\nData is not a readable file")

    if (os.path.exists(data) == False):
        print ("\nData does not exist")

    return (data, LRG_xml)


##########

# Parsing .xml file with 'xml.etree.ElementTree'

def handle_xml(data):
    """
    Parse xml files
    """

    tree = ET.parse(data)
    root = tree.getroot()
    up_anno = tree.getroot()[1]
    return(root, up_anno)


def get_gen_data(data):
    """
    Extract gene name (HGVS nomenclature) and tag strand as forward(+) or reverse(-)
    """
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


def get_background(root):
    """
    Get gene background of the LRG file, including the schema, LRG ID, HGVS ID,
    source of sequence, transcrip, coordinate system, their start, end and strand
    coordinates 
    """

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
    """
    Get information about the Genome Reference Consortium human (GRCh) build
    """

    build_lst = up_anno[1].findall('mapping')
    print('\nNumber of GRCh builds: ', len(build_lst))

    build_data = []

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
                   lrg_size = (1+(lrg_end - lrg_start))
                   lrg_size_37 = lrg_size
                else:
                   lrg_size = (1+(lrg_start - lrg_end))
                   lrg_size_37 = lrg_size

#            print('\n' + build, NC_trans, gstart, gend, lrg_start, lrg_end, lrg_size)

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
                   lrg_size = (1+(lrg_end - lrg_start))
                   lrg_size_38 = lrg_size
                else:
                   lrg_size = (1+(lrg_start - lrg_end))
                   lrg_size_38 = lrg_size

#            print('\n' + build, NC_trans, gstart, gend, lrg_start, lrg_end, lrg_size)

        # if any build other than 37 or 38 is present, this script will need to be modified
        else:
            build = annotation.get('coord_system')
            chro = annotation.get('other_name')
            NC_trans = annotation.get('other_id')
            gstart = annotation.get('other_start')
            gend = annotation.get('other_end')
            lrg_start = int(lrg.get('lrg_start'))
            lrg_end = int(lrg.get('lrg_end'))

        # create list of variables to be exported to .txt file called build_data
        build_data.append([build, chro, NC_trans, str(gstart), str(gend), str(lrg_start), str(lrg_end), str(lrg_size)])

    # set data to be exported into columns
    for group in build_data:
        pass

    return(build_data, build, chro, NC_trans, gstart, gend, lrg_start, lrg_end, lrg_size)

#    return (build, chro, NC_trans, gstart, gend, lrg_start, lrg_end, lrg_size)

    if lrg_size_38 == lrg_size_37:
        print("\nLRG size is the same between GRCh37 and GRCh38")
    else:
        print("\nWARNING: LRG sizes differ between GRCh37 and GRCh38")

#def is_int(gstart, gend):
#    assert str(gstart).isdigit()
#    assert str(gend).isdigit()
#    return(gstart_int,gend_int)


def get_exon_data(data, gstart, gend, chro, str_dir):
    """
    Get information about the exon for the different transcripts, including number
    of exons, exons coordinates in the LRG system regarding the cdna, transcript and protein
    """
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
                    g_start_ex = str(int(start_ex) + int(gstart))
                    g_end_ex = str(int(end_ex) + int(gend))

                else:
                    print ("\nProblem extracting exon information")

            # Create list of coordinates
            list4bed.append([chro, g_start_ex, g_end_ex, str_dir, str(trans_number)])
            list_all_coord.append([str(trans_number), ex_num, start_ex, end_ex, start_ex_tr, end_ex_tr, start_ex_pt,end_ex_pt])

    # Preparing lists to be printed in columns
    for group in list_all_coord:
        pass

    return (list_all_coord, list4bed)

# Parse xml and retrieve data for all sequence differences between build 37 and 38

def diff_data(data):
    """
    Get differences between annotations
    """
    diff_data = []

    for diff in root.findall('./updatable_annotation/annotation_set[@type="lrg"]/mapping[@type="main_assembly"]/mapping_span/diff'):
        diff_type = diff.get('type')
        diff_lrg_start = diff.get('lrg_start')
        diff_lrg_end = diff.get('lrg_end')
        diff_gen_start = diff.get('other_start')
        diff_gen_end = diff.get('other_end')
        lrg_seq_base = diff.get('lrg_sequence')
        other_seq_base = diff.get('other_sequence')

        # create list of differences between builds data
        diff_data.append([diff_type,diff_lrg_start,diff_lrg_end,diff_gen_start,diff_gen_end,lrg_seq_base,other_seq_base])

    # Data to be returned in columns in .csv file
    for group in diff_data:
        pass

    return(diff_data)


def coord2file(list_all_coord, list4bed):
    """
    Creates .csv, a comma separated text file with exon, transcripts, protein coordinates and a bed file
    """

    # Open file in read/write format. If one doesn't exist, it  will create a new one
    db = open(opath + lrg_id + "_" + enter_gene + "_coord.txt","w")
    db_csv = open(opath + lrg_id + "_" + enter_gene + "_coord.csv","w")
    bed = open (opath + lrg_id + "_" + enter_gene + "_bed", "w")

    # Add headers to files
    headings = ["trans_vers","exon", "exon_start", "exon_end", "trans_start", "trans_end", "prot_start", "prot_end"]
    bed_headings = ["chrom", "chromStart", "chromEnd", "strand", "transcript" ]

    # Writing tab separated text file for all coords
    db.write("\t".join(headings) + "\n") # writing headings
    for group in list_all_coord:
        db.write ("\t".join(group) + "\n") # writing coordinates

    # Writing comma-seperated values (.csv) file for all coords
    db_csv.write(",".join(headings) + "\n") # writing headings
    for group in list_all_coord:
        db_csv.write (",".join(group) + "\n") # writing coordinates

    # Writing tab-delimited .bed file
    bed.write("\t".join(bed_headings) + "\n") # writing headings
    for group in list4bed:
        bed.write ("\t".join(group) + "\n")
    
    return

def diff2file (build_data, diff_data):
    """
    Creates a a comma and tab separated file describing the build differences
    """
    
    # Open files
    db_build = open(opath + lrg_id + "_" + enter_gene + "_build.txt","w")
    db_diff = open(opath + lrg_id + "_" + enter_gene + "_diff.csv", "w")
    
    # Create headings
    build_headings = ["build", "chr","NC_trans","geno_start_coord", "geno_end_coord","LRG_start", "LRG_end","LRG_size"]
    diff_headings = ["Diff_type","LRG_start","LRG_end","Geno_coord_start","Geno_coord_end","LRG_seq_base","LRG_other_base"]

    # Writing a comma-separated values (.csv) 
    db_build.write("\t".join(build_headings) + "\n")
    for group in build_data:
        db_build.write("\t".join(group) + "\n")

    # Writing a comma-separated values (.csv) 
    db_diff.write(",".join(diff_headings) + "\n")
    for group in diff_data:
        db_diff.write(",".join(group) + "\n")
    
    # Close files
    db_build.close()
    db_diff.close()
    
    return
    
    # Inform where to find output files
    print ("\nOutput files have been saved in /Outputs!")
    
def disclaimer():
    print ("""\nPlease cite this software as: Gomez-Paramio, I. and Fryer, V. (2016), 'lrgext', Software, 
    Faculty of Medicine and Human Sciences, The University of Manchester.' or successor 
    references as defined by the authors.\n""")
    return

def final_tests():
    """
    Tests run at the end of the program
    """

    """ Checking for xml format """
    if schema != "1.9":
        print ("""\nlrgext supports LRG xmls built using schema 1.9. 
        Please be aware of possible data incongruencies""")

    # Check the strand direction and warn if in reverse.
    # Use LRG_571 for a forward strand example
    if str_dir == "-":
        print("\nN.B. This LRG is on the REVERSE strand")
    elif str_dir == '+':
        print("\nN.B. This LRG is on the FORWARD strand")
    return


#### MAIN information ####

(data, LRG_xml)=initial_tests()
create_repository_file()
#(gstart_int,gend_int) = is_int(gstart,gend)
(root, up_anno) = handle_xml(data)  

(schema, lrg_id,  hgnc_id, seq_source, transcript, cs, start_cs, end_cs, strand_cs) = get_background(root) # 3
(gene, str_dir) = get_gen_data(data) 

(build_data, build, chro, NC_trans, gstart, gend, lrg_start, lrg_end, lrg_size) = get_build_info(up_anno) # 5
(diff_data) = diff_data(data)

(list_all_coord, list4bed) = get_exon_data(data, gstart, gend, chro, str_dir) 
coord2file(list_all_coord, list4bed) 
diff2file(build_data, diff_data)

final_tests()
disclaimer() 

### End of Main ###
