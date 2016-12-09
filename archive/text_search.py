import sys, os, xml.etree.ElementTree as ET, csv

path = './LRGs/'

# Usage: python text_search.py gene (e.g. COL1A1. Can be lower or upper case or combination of both)

enter_gene = sys.argv[1].upper()

with open('gene_lrg_lst.csv', 'w') as f:
    csvfileWriter = csv.writer(f)

    for filename in os.listdir(path):
        if not filename.endswith('.xml'): continue
        fullname = os.path.join(path, filename)
        tree = ET.parse(fullname)
        root = tree.getroot()
        gene = root.find('updatable_annotation/annotation_set/lrg_locus').text

        csvfileWriter.writerow([gene, filename])
    f.close()

# enter_gene = input("Enter gene name: ").upper()

with open('gene_lrg_lst.csv', 'r') as f:
    csvfileReader = csv.reader(f)
    for row in csvfileReader:
        found = 'No'
        if enter_gene in row[0]:
            found = 'Yes'
            print(row[1])
            break
    if found == 'No':
        print("No LRG.xml can be found for gene entered.\nCheck gene has an associated LRG at www.lrg-sequence.org")
