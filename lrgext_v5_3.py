# -*- coding: utf-8 -*-
"""
@authors: Idoia Gomez-Paramio and Verity Fryer 2016

"""

import xml.etree.ElementTree as ET
import os.path   
# import sys
"""
Usage: 
lrgext extract build, gene, transcript, exon information from a file in LRG format 
producing different output files:
    - cvs file
    - tab separated txt file
    - bed file

"""    
    

"""

https://www.tutorialspoint.com/python/python_command_line_arguments.htm
getopt.getopt(LRG, -l:, [long_options])
"""
LRG = 'LRG_62'
path = './LRGs/'
data = path + LRG + '.xml'


# ask user to input LRG name
# filename = input("Enter LRG name: ")
# change to  tree = ET.parse(filename + '.xml') once program is ready

def handle_xml(data):
    """Use the 'xml.etree.ElementTreee' to extract inoformation from xml files"""
    
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
    and genomic start and end. It will provide "N/A", when protein coordinates are not
    available"""
    
    for annotation in up_anno[1].findall('mapping'):
        coord = annotation.get('coord_system')
        chro = annotation.get('other_name')
        NC_trans = annotation.get('other_id')
        gstart = annotation.get('other_start')
        gend = annotation.get('other_end')
        print (coord, chro, NC_trans, gstart, gend)
        return (coord, chro, NC_trans, gstart, gend)
        
        ####VF CODE###
        # Not sure if the differences between builds could be included here or 
        # in a different function.
        
        ## Also to capture if difference in sequences occur in the intro/exon
        
            
def get_gen_data(data):
    """Extract gene name (HGVN nomenclature) and tag strand as 
    forward(+) or reverse(-) """
    
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
    in the LRG system regarding the cdna, transcript and protein """
    
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
                    print ("problem when extrating exon info")
                    
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
    
def disclaimer():
    print ( 
    """Please cite this software as: 'Gomez-Paramio, I. and Fryer, V. (2016), 'lrgext', Software, 
    Faculty of Medicine and Human Sciences, The University of Manchester.' or successor 
    references as defined by the authors.\n""")
    
    return
    
##### TESTING #####

def initial_tests():
    """ Run initial tests to check the software and file before execution"""
    
    if (os.path.isfile(data) == False) :
        print ("Data is not a readable file")
    
    if (os.path.exists(data) == False):
        print ("Data does not exits")
    
    #### VF's code ####
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
    #### IDOIA CODE ####
    
    return

def builds():
    
    #### VF's code
    # if transcript number differ between builds 37 and 38:
    # print("Do transcripts match between builds? ", True/False)
    
    return
    
### End of testing ###   
   
#### MAIN information ####

initial_tests()   # 1
(root, up_anno) = handle_xml(data)  # 2
(schema, lrg_id,  hgnc_id, seq_source, transcript, cs, start_cs, end_cs, strand_cs) = get_background(root) # 3
(coord, chro, NC_trans, gstart, gend) = get_build_info(up_anno) # 4
(gene, str_dir) = get_gen_data(data) # 5
(list_all_coord, list4bed)= get_exon_data(data, gstart, gend, chro, str_dir) # 6
output2file(list_all_coord, list4bed) # 7
disclaimer() # 8
final_tests() # 9 


### End of Main ###

""" THIS SECTION CAN BE DELETED AT THE END

Main program running through different steps:
    1. Initial tests
    2. Handle xml structure
    3. Get background information
    4. Get information about the different builds
    5. Get information about the gene and strand
    6. Get information about the exon
    7. Save exon, transcript and protein coordinates to file
    8. Run of final tests


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
- Bug resolved for LRG lacking first prot. coordinates in exon 1 
    (e.g. LRG_292 and LRG 62). Comments revised, added and updated v.5.3


...... TO BE DONE
1. Compare builds indicating if seq. changes occur in introns/exons: VF
2. Add test to check for strand: IGP 
3. Add test to check if file is in xml format: VF
4. Add test regarding builds: VF
5. Any other tests suggested (See section Control in Readme file):VF


...... OUTLOOK (time allowed)
1. Providing comand line input to access directly the files/webpage
2. Adding new test using 'assert'
3. Automatic creation of output file names


"""
