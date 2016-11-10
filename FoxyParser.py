# Python 3.5
import xml.etree.ElementTree as ET

def xml_checker(file_name):
    '''Runs a check to confirm a .xml file has been entered'''
    assert file_name.endswith('.xml'), ' wrong file type entered'
    return(check_file(file_name))
    
def check_file(file_name):
    ''' check the input file is an LRG and return root '''    
    tree = ET.parse(file_name)
    root = tree.getroot()
    assert root.tag == 'lrg', 'Input file must be an LRG'
    return(root)

print(xml_checker('LRG_62_dud.xml'))
#print(xml_checker('wrong.txt'))