import urllib.request, urllib.error
import xmltodict

def homepage(request):
    file = urllib2.urlopen('view-source:http://ftp.ebi.ac.uk/pub/databases/lrgex/LRG_292.xml')
    data = file.read()
    file.close()

    data = xmltodict.parse(data)
    return render_to_response('my_template.html', {'data': data})


#import xml.etree.ElementTree as ET
#tree = ET.parse('view-source:http://ftp.ebi.ac.uk/pub/databases/lrgex/LRG_292.xml')
#root = tree.getroot()

#for child in root:
#    print(child.tag, child.attrib)
