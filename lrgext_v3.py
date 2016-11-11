import xml.etree.ElementTree as ET
import os.path   

data = 'LRG_1.xml'


# ask user to input LRG name
# filename = input("Enter LRG name: ")

# check file is in xml format. If not, return error message "Not an xml file"

# add try, except to close program if no LRG exists

# change to  tree = ET.parse(filename + '.xml') once program is ready



def get_structure(data):
    tree = ET.parse(data)
    root = tree.getroot()
    fix_anno = tree.getroot()[0]
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
    """Get build information """
    
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
    
    path_strand = "('./updatable_annotation/annotation_set')[2]"
    gene = root.find('updatable_annotation/annotation_set/lrg_locus').text
    print('Gene: ', gene)

    for annotation in root.findall('./updatable_annotation/annotation_set[@type="lrg"]/mapping[@type="main_assembly"]/mapping_span'):
        str_dir = annotation.attrib['strand']
        print (str_dir)
            
    return gene



def get_exon_data(data):
    """Get information about the exon, including number of exons, exons coordinates 
    in the LRG system regarding the cdna, transcript and protein. Information is included 
    in 3 different lists"""
    
    exon_lst = root.findall('./fixed_annotation/transcript/exon')
    print('Exon count: ', len(exon_lst))
    list_ex, list_ex_tr,list_ex_pt, list_ex  = [],[],[],[]
    
    
    for exons in root.findall('./fixed_annotation/transcript/exon'):
        ex_num = exons.get('label')
        list_ex.append(ex_num) 
        
        for exons in root.findall('./fixed_annotation/transcript/exon/coordinates'):
            
            if (exons.get('coord_system').find("t1")!=-1):
                start_ex_tr = exons.get('start')
                end_ex_tr = exons.get('end')
                list_ex_tr.append([ex_num, start_ex_tr, end_ex_tr])
            
            elif (exons.get('coord_system').find("p1")!=-1):
                start_ex_pt = exons.get('start')
                end_ex_pt = exons.get('end')
                list_ex_pt.append([ex_num, start_ex_pt, end_ex_pt])
                
            else:
                start_ex = exons.get('start')
                end_ex = exons.get('end')
                list_ex.append([ex_num, start_ex, end_ex])
                
    return (list_ex, list_ex_tr, list_ex_pt)
     
def printing():
    """ Printing of functions' outputs """
    
    print ("\nExons coordinates: ", list_ex)
    print ("\nTranscript coordinates: ", list_ex_tr)
    print ("\nProtein coordinates: ", list_ex_pt)
    
    #db = open("db.txt","r+")
    #db.write("Holaaaaa", list_ex)
    #db.close()
    
    return      
db = open("db.txt","r+")
db.write("Holaaaaa")
db.close()

    
# if transcript number differ between builds 37 and 38:
# print("Do transcripts match between builds? ", True/False)

# ouput all to .csv file or BED file
# tree.write('output.txt')


##### TESTING #####
""" Section dedicated to testing the running of the programm. It is subdivided into
initial and secondary tests"""

def initial_tests():
    """ Tests run at the beggining of the program"""
    
    if (os.path.isfile(data) == False) :
        print ("Data is not a readable file")
    
    if (os.path.exists(data) == False):
        print ("Data does not exits")
    return
    
def final_tests():  
    """Tests run at the end of the program"""
    
    if schema != "1.9":  # Checking for xml format
        print ("""lrgext supports LRG xmls built using schema 1.9. Please" 
        be aware of data incongruencies""")
    return
    
def strand_dir():
    
    
    return
### End of testing ###   
    
#### MAIN ####
""" Main program running through different steps:
    1. Initial tests
    2. Getting background information
    3. Get information about the different builds
    4. Get information about the strand
    5. Get information about the exon
    6. Print in the standard output
    7. Running of the final tests
"""

initial_tests()   # 1
(root, up_anno) = get_structure(data)  # 2
(schema, lrg_id,  hgnc_id, seq_source, transcript, cs, start_cs, end_cs, strand_cs) = get_background(root) # 3
(coord, chro, NC_trans, gstart, gend) = get_build_info(up_anno) # 4
(gene) = get_up_anno(data) # 5
(list_ex, list_ex_tr, list_ex_pt)= get_exon_data(data) # 6
printing() # 7
final_tests()  # 8

### End of Main ###

""" to be done
1. Compare builds
2. Output, generate files
3. Deal with issues, as having two transcripts
"""

""" Outlook 
1. Providing comand line input to access directly the webpage
2. Adding new test using 'assert'
"""
