# this could be used to extract info from a webpage. Not necessary because labs will usually have the info saved

#from urllib.request import urlopen

#LRG = "LRG_1"
#page = "http://ftp.ebi.ac.uk/pub/databases/lrgex/" + LRG + ".xml"

#URL_xml = urlopen(page)
#xml = URL_xml.read()
#print (xml)

import xml.etree.ElementTree as ET

# add try, except to close program if no LRG exists

tree = ET.parse('LRG_1.xml')
root = tree.getroot()
fix_anno = tree.getroot()[0]
up_anno = tree.getroot()[1]

# Print LRG_id, NG, dna source etc.
print ("LRG id = ", root[0][0].text)
print ("Sequence source = ", root[0][2].text)
print ("Organsim = ", root[0][3].text)
print ("Source of data = ", root[0][4][0].text)
print ("Mol_type = ", root[0][5].text)
print ("LRG creation date = ", root[0][6].text)

#for annotation in up_anno.findall('annotation_set'):
for annotation in up_anno[1].findall('mapping'):
    #GHC = annotation.find('rank').text
    #name = annotation.get('type')
    coord = annotation.get('coord_system')
    chro = annotation.get('other_name')
    NC_trans = annotation.get('other_id')
    gstart = annotation.get('other_start')
    gend = annotation.get('other_end')
    print ("Build",coord, chro, NC_trans, gstart, gend, )

# if transcript number differ between builds 37 and 38:
# print("Do transcripts match between builds? ", True/False)

# Find strand of LRG (forward 1, reverse -1)
for annotation in up_anno[1].findall('mapping_span'):
    #GHC = annotation.find('rank').text
    #name = annotation.get('type')
    strand = annotation.get('strand')
    print ("Strand = ", strand)

# Print gene name
#? print(tree.findall[0]('symbol'))
#? gene = gene_name.get('symbol')
#? gene = tree.getroot()[1][3][16][0][0]
#? print ("Gene = ", gene.attrib)

for exon in fix_anno.iter('exon'):
    print (exon.attrib)

# ouput all to .csv file or BED file
