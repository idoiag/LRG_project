
# this could be used to extract info from a webpage. Not necessary because labs will usually have the info saved


import xml.etree.ElementTree as ET
import os.path   ### NEW

data = 'LRG_517.xml'

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
    chro37 = tree.getroot()[1][1][2]
    chro38 = tree.getroot()[1][1][3]
    return(root, up_anno)

def get_background(root):
    
    for fixed in root.findall("fixed_annotation"):
        lrg_id = fixed.find('id').text
        hgnc_id = fixed.find ('hgnc_id').text
        seq_source = fixed.find ('sequence_source').text

        for transcript in root.findall("fixed_annotation/transcript"):
            transcript = transcript.get('name')

            for coordinates in root.findall("fixed_annotation/transcript/coordinates"):
                cs = coordinates.get('coord_system')
                start_cs = coordinates.get('start')
                end_cs = coordinates.get('end')
                strand_cs = coordinates.get('strand')

        print ( lrg_id,  hgnc_id, seq_source, transcript, cs, start_cs, end_cs, strand_cs)
        return ( lrg_id,  hgnc_id, seq_source, transcript, cs, start_cs, end_cs, strand_cs)

def get_up_anno(data):
    gene = tree.find('updatable_annotation/annotation_set/lrg_locus').text
    print('Gene: ', gene)
    return gene

def get_background(root):
    for lrg in root.findall ("."):  
        schema = lrg.get('schema_version')  
        
    for fixed in root.findall("./fixed_annotation"):
        lrg_id = fixed.find('id').text
        hgnc_id = fixed.find ('hgnc_id').text
        seq_source = fixed.find ('sequence_source').text

        for transcript in root.findall("./fixed_annotation/transcript"):
            transcript = transcript.get('name')

            for coordinates in root.findall("./fixed_annotation/transcript/coordinates"):
                cs = coordinates.get('coord_system')
                start_cs = coordinates.get('start')
                end_cs = coordinates.get('end')
                strand_cs = coordinates.get('strand')
        
        print (schema) 
        print ( lrg_id,  hgnc_id, seq_source, transcript, cs, start_cs, end_cs, strand_cs)
        return ( schema, lrg_id,  hgnc_id, seq_source, transcript, cs, start_cs, end_cs, strand_cs)


def get_build_info(up_anno):
    for annotation in up_anno[1].findall('mapping'):
        coord = annotation.get('coord_system')
        chro = annotation.get('other_name')
        NC_trans = annotation.get('other_id')
        gstart = annotation.get('other_start')
        gend = annotation.get('other_end')
        print (coord, chro, NC_trans, gstart, gend)
        return (coord, chro, NC_trans, gstart, gend)

"""
def exon_lst(filename):
    exon_lst = lrg_xml.findall('fixed_annotation/transcript/exon')
    print('Exon count: ', len(exon_lst))
    for exons in exon_lst:
        #print('Exon number: ', exons.get("label"))
        exon_number = exons.get("label")
        if int(exon_number) > 0:
            for coordinates in exons:
                if coordinates.get('coord_system') == lrg_id:
                    #print('LRG start: ', coordinates.get("start"))
                    coord_start = coordinates.get("start")
                    #print('LRG end: ', coordinates.get("end"))
                    coord_end = coordinates.get("end")
                    #print('Strand: ', coordinates.get("strand"))
                    strand = coordinates.get("strand")
        print(exon_number, coord_start, coord_end, strand)

"""

# if transcript number differ between builds 37 and 38:
# print("Do transcripts match between builds? ", True/False)

# ouput all to .csv file or BED file
# tree.write('output.txt')


#### TESTING START ########
def initial_tests():
    if (os.path.isfile(data) == False) :
        print ("Data is not a readable file")
    
    if (os.path.exists(data) == False):
        print ("Data does not exits")
    
def second_tests():  
    if schema != "1.9":  # Checking for xml format
        print ("""lrgext supports LRG xmls built using schema 1.9. Please" 
        be aware of data incongruencies""")
    return

##### END OF TESTING ########    
    
#### MAIN ####
initial_tests()   
(root, up_anno) = get_structure(data)  
(schema, lrg_id,  hgnc_id, seq_source, transcript, cs, start_cs, end_cs, strand_cs) = get_background(root)
(coord, chro, NC_trans, gstart, gend) = get_build_info(up_anno)
(gene) = get_up_anno(data)
second_tests()  

#### Not working

#(exon_number, coord_start, coord_end, strand) = exon_lst(filename)
