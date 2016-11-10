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
    return(check_public(root))

def check_public(root):
    ''' Checks that the LRG file is a public file. for pending files issues a warning regarding completeness'''
    if root.findall("./fixed_annotation/.")[0].tag == 'source':
        return(root)
    else:
        print ('Warning! this is a pending LRG file and may be subject to modification')
        return(root)

print(xml_checker('LRG_62.xml'))
#print(xml_checker('wrong.txt'))