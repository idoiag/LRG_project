# -*- coding: utf-8 -*-
"""
@authors: Idoia Gomez-Paramio and Verity Fryer 2016

"""

import xml.etree.ElementTree as ET
import os.path   
# import sys

#getopt.getopt(LRG, -l:, [long_options])

"""
https://www.tutorialspoint.com/python/python_command_line_arguments.htm
"""
LRG = 'LRG_517'
path = './LRGs/'
data = path + LRG + '.xml'


# ask user to input LRG name
# filename = input("Enter LRG name: ")

# change to  tree = ET.parse(filename + '.xml') once program is ready

def get_structure(data):
    tree = ET.parse(data)
    root = tree.getroot()
    up_anno = tree.getroot()[1]
    return(root, up_anno)

def get_background(root):
    """Get background information about the gene """
    
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
        
        print (schema) 
        print ( lrg_id,  hgnc_id, seq_source, transcript, cs, start_cs, end_cs, strand_cs)
        return ( schema, lrg_id,  hgnc_id, seq_source, transcript, cs, start_cs, end_cs, strand_cs)


def get_build_info(up_anno):
    """Get build information, including coordinates, chromosome, transcript,
    and genomic start and end"""
    
    for annotation in up_anno[1].findall('mapping'):
        coord = annotation.get('coord_system')
        chro = annotation.get('other_name')
        NC_trans = annotation.get('other_id')
        gstart = annotation.get('other_start')
        gend = annotation.get('other_end')
        print (coord, chro, NC_trans, gstart, gend)
        return (coord, chro, NC_trans, gstart, gend)
            
def get_up_anno(data):
    """Get information about the strand information """
    
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
    
def get_exon_data(data, gstart, gend, chro, str_dir):
    """Get information about the exon for the different transcripts, including number of exons, exons coordinates 
    in the LRG system regarding the cdna, transcript and protein. """
    
    trans_number = 0
    list_all_coord, list4bed = [], []
    
    # Loop to extract information when there is more than 1 transcript 
    # (e.g.LRG_214 has two transcripts
    for transcripts in root.findall('./fixed_annotation/transcript'):
        trans_number += 1
     
        exon_lst = root.findall('./fixed_annotation/transcript/exon')
        print( "Transcript number: ", trans_number, 'Exon count: ', len(exon_lst))
    
        for exons in transcripts.findall('exon'):
            ex_num = exons.get('label')
            
            # Extracting transcript, protein and exon coordinates (in this order)
            for coord in exons.findall('coordinates'):
                
                if (coord.get('coord_system').find("t")!=-1):
                    start_ex_tr = coord.get('start')
                    end_ex_tr = coord.get('end')
                
                elif (coord.get('coord_system').find("p")!=-1):
                    start_ex_pt = coord.get('start')
                    end_ex_pt = coord.get('end')
                    
                else:
                    start_ex = coord.get('start')
                    end_ex = coord.get('end')
                    g_start_ex = start_ex + gstart
                    g_end_ex = end_ex + gend
            
            # Create list of coordinates
            list4bed.append([chro, g_start_ex, g_end_ex, str_dir, str(trans_number)])
            list_all_coord.append([str(trans_number), ex_num, start_ex, end_ex, start_ex_tr, end_ex_tr, start_ex_pt,end_ex_pt]  )
    
    # Preparing lists to be print in columns
    for group in list_all_coord: 
        print ("\t".join(group) + "\n")
                
    return (list_all_coord, list4bed)
     
def output2file(list_all_coord, list4bed):
    """ Creating csv, a tab separate txt file with exon, transcripts and 
    protein coordinates and a bed file"""

    # Opening files to read and write. If not existing, it will create a new one
    db = open("./Outputs/LRG_coord.txt","w")
    db_csv = open("./Outputs/LRG_coord.csv","w")
    bed = open ("./Outputs/LRG_bed", "w")
    
    # Add headings to files
    headings = ["transcript","exon", "ex_start", "ex_end", "tr_start", "tr_end", "pt_start", "pt_end"]
    bed_headings = ["chr", "start", "end", "strand", "transcript" ]
    
    #Writting tab separated text file
    db.write("\t".join(headings) + "\n") # writting headings
    for group in list_all_coord:   
        db.write ("\t".join(group) + "\n") # writting coordinates
    
    #Writting csv file
    db_csv.write(",".join(headings) + "\n") # writting headings
    for group in list_all_coord:   
        db_csv.write (",".join(group) + "\n") # writting coordinates
    
    #Writting bed file
    bed.write("\t".join(bed_headings) + "\n") # writting headings
    for group in list4bed:   
        bed.write ("\t".join(group) + "\n")
    
    # Closing files
    db.close()
    db_csv.close()
    bed.close()
    
    return      
    
def creating_bed():
    
    
    return    
    
    
def disclaimer():
    print ( 
    """Please cite this software as: 'Gomez-Paramio, I. and Fryer, V. (2016), 'lrgext', Software, 
    Faculty of Medicine and Human Sciences, The University of Manchester.' or successor 
    references as defined by the authors.\n""")
    
    return
    
##### TESTING #####
""" Section dedicated to testing the running of the programm. It is subdivided into
initial and secondary tests"""

def initial_tests():
    """ Tests run at the beggining of the program"""
    
    if (os.path.isfile(data) == False) :
        print ("Data is not a readable file")
    
    if (os.path.exists(data) == False):
        print ("Data does not exits")
     
    # check file is in xml format. If not, return error message "Not an xml file"
    # add try, except to close program if no LRG exists
        
    return
    
    
def final_tests():  
    """Tests run at the end of the program"""
    
    if schema != "1.9":  # Checking for xml format
        print ("""lrgext supports LRG xmls built using schema 1.9. Please" 
        be aware of data incongruencies""")
    return
    
def strand_dir():
    """Check the strand direction and warn if in reverse.
    Use LRG_571 for a forward strad example  """
    
    return

def builds():
    # if transcript number differ between builds 37 and 38:
    # print("Do transcripts match between builds? ", True/False)
    return
    
### End of testing ###   
   
#### MAIN ####
""" Main program running through different steps:
    1. Initial tests
    2. Getting background information
    3. Get information about the different builds
    4. Get information about the strand
    5. Get information about the exon
    6. Save exon, transcrip and protein coordinates to file
    7. Running of the final tests
"""

initial_tests()   # 1
(root, up_anno) = get_structure(data)  # 2
(schema, lrg_id,  hgnc_id, seq_source, transcript, cs, start_cs, end_cs, strand_cs) = get_background(root) # 3
(coord, chro, NC_trans, gstart, gend) = get_build_info(up_anno) # 4
(gene, str_dir) = get_up_anno(data) # 5
(list_all_coord, list4bed)= get_exon_data(data, gstart, gend, chro, str_dir) # 6
output2file(list_all_coord, list4bed) # 7
disclaimer() # 8
final_tests() # 9 


### End of Main ###
""" 
.......DONE: Versions
- Extracting background information about the gene: v1
- Making the code modular: v2
- Adding examples of initial and final test: v2.1
- Extracting exon, transcript and protein coordinates: v3
- Output: it creates a csv and tab separate text file: v3.1
- Dealing with more than one transcript v.4
- Adding disclaimer v.4.1
- Creation of a bed file v.5
- Organisation of input and outputs into folders v.5.1
- Splitting of input file name. Code clean, comments added v.5.2


...... TO BE DONE
1. Compare builds
2. Adding further tests
3. Resolve issues with LRG_62 and LRG_292


...... OUTLOOK (time allowed)
1. Providing comand line input to access directly the files/webpage
2. Adding new test using 'assert'
3. Automatic creation of output file names


"""
