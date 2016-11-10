# Python 3.5
import xml.etree.ElementTree as ET


def check_file(file_name):
    ''' check the input file is an LRG '''
    tree = ET.parse(file_name)
    root = tree.getroot()
    
    assert = root.tag == 'LRG', 'Input file must be an LRG - do it again doofus'
    
    return root

check_file('LRG_62.xml')