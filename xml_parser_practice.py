import xml.etree.ElementTree as ET
lrg_xml = ET.parse('LRG_292.xml')
#root = tree.getroot()

lrg_id = lrg_xml.find('fixed_annotation/id').text
print('LRG id: ', lrg_id)

def get_gene_name(filename):
    gene = lrg_xml.find('updatable_annotation/annotation_set/lrg_locus').text
    print('Gene: ', gene)
    return gene

get_gene_name('LRG_292.xml')

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


#    output.write(exon_number,'\t',coord_start,'\t',coord_end,'\t',strand,'\n')

#import urllib.request, urllib.error
#import xmltodict

#def homepage(request):
#    file = urllib2.urlopen('view-source:http://ftp.ebi.ac.uk/pub/databases/lrgex/LRG_292.xml')
#    data = file.read()
#    file.close()

#    data = xmltodict.parse(data)
#    return render_to_response('my_template.html', {'data': data})


#import xml.etree.ElementTree as ET
#tree = ET.parse('view-source:http://ftp.ebi.ac.uk/pub/databases/lrgex/LRG_292.xml')
#root = tree.getroot()

#for child in root:
#    print(child.tag, child.attrib)
