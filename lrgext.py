# -*- coding: utf-8 -*-

import xml.etree.ElementTree as ET

def get_structure(data):
    tree = ET.parse(data)
    root = tree.getroot()
    fix_anno = tree.getroot()[0]
    up_anno = tree.getroot()[1]
    chro37 = tree.getroot()[1][1][2]
    chro38 = tree.getroot()[1][1][3]
    return(root, up_anno)

def get_background(root):
    
    for fixed in root.findall("./fixed_annotation"):
        lrg_id = fixed.find('id').text
        hgnc_id = fixed.find ('hgnc_id').text
        seq_source = fixed.find ('sequence_source').text
        
        for transcript in root.findall("./fixed_annotation/transcript"):
            transcript = transcript.get('name')
            
            #### Variables still need to be included in the return/print statement
            
            for coordinates in root.findall("./fixed_annotation/transcript/coordinates"):
                cs = coordinates.get('coord_system')
                start_cs = coordinates.get('start')
                end_cs = coordinates.get('end')
                strand_cs = coordinates.get('strand')
                
            for coordinates in root.findall("./fixed_annotation/transcript/coding_region/coordinates"):
                cs_cr = coordinates.get('coord_system')
                start_cs_cr = coordinates.get('start')
                end_cs_cr = coordinates.get('end')
                strand_cs_cr = coordinates.get('strand')
                
            for coor_exon in root.findall("./fixed_annotation/transcript/exon"):
                exon_number = coordinates.get('label')
                
                for coor_exon in root.findall("./fixed_annotation/transcript/exon/coordinates"):
                    exonC_start = coor_exon.get('start')
                    exonC_end = coor_exon.get('end')
                    
            ###### End of variables ########
                
        print ( lrg_id,  hgnc_id, seq_source, transcript, cs, start_cs, end_cs, strand_cs)
        return ( lrg_id,  hgnc_id, seq_source, transcript, cs, start_cs, end_cs, strand_cs)
        #return

def get_build_info(up_anno):        
    for annotation in up_anno[1].findall('mapping'):
        coord = annotation.get('coord_system')
        chro = annotation.get('other_name')
        NC_trans = annotation.get('other_id')
        gstart = annotation.get('other_start')
        gend = annotation.get('other_end')
        print (coord, chro, NC_trans, gstart, gend, )
        return (coord, chro, NC_trans, gstart, gend)
        #return

        
#### MAIN ####
(root, up_anno) = get_structure('LRG_1.xml')
(lrg_id,  hgnc_id, seq_source, transcript, cs, start_cs, end_cs, strand_cs) = get_background(root)
(coord, chro, NC_trans, gstart, gend) = get_build_info(up_anno)
