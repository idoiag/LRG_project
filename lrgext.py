# -*- coding: utf-8 -*-
import xml.etree.ElementTree as ET, os.path, sys, csv, getopt
"""
Authors: Idoia Gomez-Paramio and Verity Fryer 2016
lrgext program extracts build, gene, transcript, exon information from a file in LRG format given a HGVS gene (case insensitive).
Data is exported to comma seperated values (.csv) and tab separated file (as .txt).
A bed file will also be created.
The script supports LRG files with more than one transcript.
Usage: python lrgext.py -g <gene>
"""

def create_repository_file(path, lrg_list):
    """
    This function will parse all .xml files in LRGs directory ("LRGs") and create a csv file.
    This file ("gene_lrg_lst.csv") lists all existing .xml files within the directory, and the gene name associated with the LRG
    """

    with open(lrg_list,'w',  newline='') as f:
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

def handle_xml(data):
    """
    Function to parse .xml file with 'xml.etree.ElementTree'
    """

    tree = ET.parse(data)
    root = tree.getroot()
    up_anno = tree.getroot()[1]
    return(root, up_anno)

def get_gen_data(data, root, lrg_id):
    """
    Extract gene name (HGVS nomenclature) and tag strand as forward(+) or reverse(-)
    """

    # Extract name of gene
    gene = root.find('updatable_annotation/annotation_set/lrg_locus').text
    #print('\nGene: ', gene)
    print ("\nLRG ID: ", lrg_id, ' Gene: ', gene)

    # Determine the strand direction of gene, forward strand or reverse strand
    for annotation in root.findall('./updatable_annotation/annotation_set[@type="lrg"]/mapping[@type="main_assembly"]/mapping_span'):
        str_dir = annotation.attrib['strand']

        # Mark forward strand as "+" and negative as "-"
        if (str_dir == "1"):
            str_dir = "+"
        else:
            str_dir = "-"

    return (gene, str_dir)

def get_background(root):
    """
    Get gene background of the LRG file, including xml schema, LRG ID, HGVS ID, source of sequence, transcript, coordinate system.
    Also capture start and end co-ordinates and strand. 
    It will return "N/A", when protein coordinates are not available
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
        #print ("\nLRG ID: " + lrg_id)
        return (schema, lrg_id, hgnc_id, seq_source, transcript, cs, start_cs, end_cs, strand_cs)

def get_build_info(up_anno):

    """
    Capture information about the Genome Reference Consortium human (GRCh) build
    """

    build_lst = up_anno[1].findall('mapping')
    print('\nNumber of GRCh builds: ', len(build_lst), "\n")

    build_data = []

    for annotation in up_anno[1].findall('mapping'):
        build = annotation.get('coord_system')
        chro = annotation.get('other_name')

        # collect info on build 37
        if build.startswith('GRCh37'):
            NC_trans = annotation.get('other_id')
            gstart = int(annotation.get('other_start'))
            gend = int(annotation.get('other_end'))
            gene_len = (gend - gstart) 
            gene_len_37 = gene_len

            # determine start and end of LRG
            for lrg in up_anno[1].findall('mapping/mapping_span'):
                lrg_start = int(lrg.get('lrg_start'))
                lrg_end = int(lrg.get('lrg_end'))
                # check if lrg end is larger than lrg start (which it should be regardless of strand direction)
                if lrg_end > lrg_start:
                   lrg_size = (lrg_end - lrg_start)
                   lrg_size_37 = lrg_size
                else:
                   lrg_size = (lrg_start - lrg_end)
                   lrg_size_37 = lrg_size

        # otherwise collect info on build 38
        elif build.startswith('GRCh38'):
            NC_trans = annotation.get('other_id')
            gstart = int(annotation.get('other_start'))
            gend = int(annotation.get('other_end'))
            gene_len = (gend - gstart) 
            gene_len_38 = gene_len

            # determine start and end of LRG
            for lrg in up_anno[1].findall('mapping/mapping_span'):
                lrg_start = int(lrg.get('lrg_start'))
                lrg_end = int(lrg.get('lrg_end'))
                # check if lrg end is larger than lrg start (which it should be regardless of strand direction)
                if lrg_end > lrg_start:
                   lrg_size = (lrg_end - lrg_start)
                   lrg_size_38 = lrg_size
                else:
                   lrg_size = (lrg_start - lrg_end)
                   lrg_size_38 = lrg_size

        # if any build other than 37 or 38 is present, this script will need to be modified to compare LRG size
        else:
            build = annotation.get('coord_system')
            chro = annotation.get('other_name')
            NC_trans = annotation.get('other_id')
            gstart = int(annotation.get('other_start'))
            gend = int(annotation.get('other_end'))
            gene_len = (gend - gstart)
            lrg_start = int(lrg.get('lrg_start'))
            lrg_end = int(lrg.get('lrg_end'))
            lrg_size = (lrg_end - lrg_start)

     # create list of variables to be exported to .txt file called build_data
        build_data.append([build, chro, NC_trans, str(gstart), str(gend), str(lrg_start), str(lrg_end), str(lrg_size), str(gene_len)])

    # set data to be exported into columns
    for group in build_data:
        pass

    return(build_data, build, chro, NC_trans, gstart, gend, lrg_start, lrg_end, lrg_size_37, lrg_size_38, lrg_size, gene_len_37, gene_len_38, gene_len)

def get_exon_data(data, gstart, gend, chro, str_dir, root):

    """
    Capture exon information for all transcripts.
    This will include number of exons and exon coordinates in the LRG system regarding the cdna, transcript and protein
    """

    trans_number = 0 # starting trans counter
    count_ex_all = 0 # starting all_exons (in all transcripts counter)
    list_all_coord, list4bed = [], []

    # Loop to extract information when there is more than 1 transcript
    for transcripts in root.findall('./fixed_annotation/transcript'):
        trans_number += 1
        count_ex_tran = 0 # to count exons per transcript

        exon_lst = root.findall('./fixed_annotation/transcript/exon')
        tot_exons = len(exon_lst)

        for exons in transcripts.findall('exon'):
            ex_num = exons.get('label')

            # Extract transcript, protein and exon coordinates (in this order)
            # Print "N/A" if protein coordinates are not available
            start_ex_pt = "N/A"
            end_ex_pt = "N/A"

            count_ex_tran +=1

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

        count_ex_all += count_ex_tran # all all exons from all transcripts
        print("Transcript: ",trans_number,  " Exons: ", str(count_ex_tran)  )

    # Prepare lists to be printed in columns
    for group in list_all_coord:
        pass

    return (list_all_coord, list4bed, tot_exons, count_ex_tran, count_ex_all)

def get_diff_data(data, root):
    """
    Parse xml and retrieve data for all sequence differences between build 37 and 38
    """
    diff_data = []

    for mapping in root.findall('./updatable_annotation/annotation_set[@type="lrg"]/mapping'):
        build = mapping.get('coord_system')
        for diff in mapping.findall('mapping_span/diff'):
            diff_type = diff.get('type')
            diff_lrg_start = diff.get('lrg_start')
            diff_lrg_end = diff.get('lrg_end')
            diff_gen_start = diff.get('other_start')
            diff_gen_end = diff.get('other_end')
            lrg_seq_base = diff.get('lrg_sequence')
            other_seq_base = diff.get('other_sequence')

            # create list of differences between builds data
            diff_data.append([build,diff_type,diff_lrg_start,diff_lrg_end,diff_gen_start,diff_gen_end,lrg_seq_base,other_seq_base])

    # Data to be returned in columns in .csv file
    for group in diff_data:
        pass

    return(diff_data)

def coord2file(opath, enter_gene, list_all_coord, list4bed, lrg_id):
    """
    Creates .csv, a comma separated text file with exon, transcripts, protein coordinates and a bed file
    """

    # Open file in read/write format. If one doesn't exist, it  will create a new one
    db = open(opath + lrg_id + "_" + enter_gene + "_coord.txt","w")
    db_csv = open(opath + lrg_id + "_" + enter_gene + "_coord.csv","w")
    bed = open (opath + lrg_id + "_" + enter_gene + ".bed", "w")

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

def diff2file (opath, enter_gene, build_data, diff_data, lrg_id):
    """
    Creates a a comma and tab separated file describing the build differences
    """

    # Open files
    db_build = open(opath + lrg_id + "_" + enter_gene + "_build.txt","w")
    db_diff = open(opath + lrg_id + "_" + enter_gene + "_diff.csv", "w")

    # Create headings
    build_headings = ["build", "chr","NC_trans","geno_start_coord", "geno_end_coord","LRG_start", "LRG_end","LRG_size", "Gene_size"]
    diff_headings = ["Build","Diff_type","LRG_start","LRG_end","Geno_coord_start","Geno_coord_end","LRG_allele","Ref_allele"]

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

    # Message for user that output files have been successfully created and where to find them
    print ("\nOutput files have been saved in /Outputs!")

def initial_tests(path, enter_gene, lrg_list):
    """
    Run initial tests to check the software and file before execution:
    1. The 'gene_lrg_lst.csv' will be searched for the gene name entered.
       If the file exist, the LRG ID will be captured. If not, a warning will be displayed.
    2. Data consistency within the file will be checked.
    """

    # Checking if gene is in list therefore if .xml file exists, capture name of LRG.xml file
    with open(lrg_list,'r') as f:
        csvfileReader = csv.reader(f)
        # Iterate through .csv file and check first column (row[0]) in each row for gene name
        # If gene name found, capture the LRG_ID from the second column (row[1)]
        for row in csvfileReader:
            found = False
            if enter_gene == row[0]:
                found = True
                LRG_xml = row[1]
                # if the gene name entered is in the first entry for a row, exit the loop
                break

        if found == False:
            # If no LRG.xml found for the gene entered, display an error message and exit the script immediately
            print("\nNo LRG.xml can be found for gene entered.\nCheck that gene has an associated LRG at www.lrg-sequence.org\n")
            exit()

    data = path + LRG_xml

    # Check data consistency within file 
    if (os.path.isfile(data) == False) :
        print ("\nData is not a readable file")

    if (os.path.exists(data) == False):
        print ("\nData does not exist")

    return (data, LRG_xml)
    
def final_tests(tot_exons, schema,  str_dir, lrg_size_37, lrg_size_38, gene_len_37, gene_len_38, count_ex_all ):
    """
    Tests run at the end of the program
    """

    # Checking for xml format
    if schema != "1.9":
        print ("""\nlrgext supports LRG xmls built using schema 1.9. 
        Please be aware of possible data incongruencies""")

    # Check the strand direction and warn if in reverse.
    if str_dir == "-":
        print("\nN.B.\t-This LRG is on the REVERSE strand")
    elif str_dir == '+':
        print("\nN.B.\t-This LRG is on the FORWARD strand")

    # check the size of the LRG is the same when comparing builds (this should ALWAYS be the same)
    if gene_len_38 != lrg_size_38:
        print("\t-GRCh38 reference sequence size differ from LRG size. Please refer to diff.csv in /Outputs folder for more information") 
    else:
        print("\t-GRCh38 reference sequence and LRG have same size") 
    
    ### Assert functions used for debugging ###
    assert (tot_exons == count_ex_all), "Problem with exon number" 
    assert (lrg_size_38 == lrg_size_37), "Warning: LRG sizes differ between builds"  

    return

def disclaimer():
    print ("""\nPlease cite this software as: Gomez-Paramio, I. and Fryer, V. (2016), 'lrgext', Software, 
Faculty of Medicine and Human Sciences, The University of Manchester.' or successor 
references as defined by the authors.\n""")
    
    return

def main(argv):
    
    ## GET OPTS ##
    enter_gene = ''
    path = ''
    opath = ''
    lrg_list = ''
    
    try:
        opts, args = getopt.getopt(argv,"si:o:g:l:",["gene=", "ifile=", 'ofile=, lrg_list='])
    except getopt.GetoptError:
        print ('lrgext.py -g <gene> -i <inputfile> -o <outputfile> -s for settings' )
        sys.exit(2)
    
    for opt, arg in opts:
        if opt == '-s':  
            print ("default settings: inputfile=./LRGs/, output=./Outputs/, LRG list =gene_lrg_lst.csv ")
            sys.exit()
        
        elif opt in ('-g', "--gene"):
            enter_gene = arg
        elif opt in ("-i", "--ifile"):
            path = arg
        elif opt in ("-o", "--ofile"):
            opath = arg
        elif opt in ("-l", "--lrg"):
            lrg_list = arg
    if enter_gene == '':
        print ("Provide name of gene -g <gene>, use -i <inputfile> -o <output> -l <LRG list> -s for settings")
    else:
        print ('Gene= ', enter_gene)  
    
    if path == '':    
        path = './LRGs/'
        print ('Default LRG_folder=', path)
    else: 
        print ('LRG_folder=', path)
    
    if opath == '':    
        opath = './Outputs/'
        print ('Default output folder=', opath)
    else: 
        print ('Output folder=', opath)
    
    if lrg_list == '':    
        lrg_list = 'gene_lrg_lst.csv'
        print ('Default LRG list=', lrg_list)
    else: 
        print ('LRG list=', lrg_list)    
        

    #lrg_list = 'gene_lrg_lst.csv'
    enter_gene=enter_gene.upper()
    
    create_repository_file(path, lrg_list)
    (data, LRG_xml)=initial_tests(path, enter_gene, lrg_list)
    (root, up_anno) = handle_xml(data)
    
    (schema, lrg_id, hgnc_id, seq_source, transcript, cs, start_cs, end_cs, strand_cs) = get_background(root)
    (gene, str_dir) = get_gen_data(data, root, lrg_id)
    
    (build_data, build, chro, NC_trans, gstart, gend, lrg_start, lrg_end, lrg_size_37, lrg_size_38, lrg_size, gene_len_37, gene_len_38, gene_len) = get_build_info(up_anno)
    (diff_data) = get_diff_data(data, root)
    
    (list_all_coord, list4bed, tot_exons, count_ex_tran, count_ex_all) = get_exon_data(data, gstart, gend, chro, str_dir, root)
    coord2file(opath, enter_gene, list_all_coord, list4bed, lrg_id)
    diff2file(opath, enter_gene, build_data, diff_data, lrg_id)
    
    final_tests(tot_exons, schema,  str_dir, lrg_size_37, lrg_size_38, gene_len_37, gene_len_38, count_ex_all)
    disclaimer()
    
if __name__ == "__main__":
    main(sys.argv[1:])


