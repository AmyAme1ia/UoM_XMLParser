    # Python 3.5
import xml.etree.ElementTree as ET
import pandas as pd

def xml_checker(file_name):
    '''Runs a check to confirm a .xml file has been entered'''
    assert file_name.endswith('.xml'), ' wrong file type entered'
    return(check_file(file_name))
    
def check_file(file_name):
    ''' check the input file is an LRG and return root '''    
    tree = ET.parse(file_name)
    root = tree.getroot()
    assert root.tag == 'lrg', 'Input file must be an LRG'
    return(get_data(root))

def check_public(root):
    ''' Checks that the LRG file is a public file. for pending files issues a warning regarding completeness'''
    if root.findall("./fixed_annotation/*")[4].tag == 'source':
        return(get_data(root))
    else:
        print ('Warning! this is a pending LRG file and may be subject to modification')
        return(get_data(root))


def get_data(root):
    ''' extract data from xml and store in a pandas dataframe (and other variables) '''
    # define an empty dataframe to accept parsed exon data relative to lrg coordinate from xml file
    df_exon_rel = pd.DataFrame(columns=['exon_no','start','end'])
    # get LRG id and genomic sequence
    lrg_id = root.findall('./fixed_annotation/id')[0].text
    genomic_sequence = root.findall('./fixed_annotation/sequence')[0].text
    # define empty lists to hold parsed exon data
    ex_label = []
    ex_start = []
    ex_end = []
    ex_strand = []
    
    # loop through exons pulling out exon number, start and stop position
    for item in root.findall('./fixed_annotation/transcript/exon'):
        ex_label.append(item.attrib['label'])
        ex_start.append(item[0].attrib['start'])
        ex_end.append(item[0].attrib['end'])
        #ex_strand.append(item[0].attrib['strand'])
    
    # enter data from lists into pandas dataframe
    for i in range(len(ex_label)):
        df_exon_rel.loc[df_exon_rel.shape[0]] = [ex_label[i],ex_start[i],ex_end[i]]

    #df_exon_rel['seq'] = genomic_sequence[df_exon_rel.start:df_exon_rel.end]
    
    
    # check that coordinates are in correct spatial orientation
    for j in range(len(df_exon_rel['start'])):
        assert int(df_exon_rel['end'].loc[j]) > int(df_exon_rel['start'].loc[j]), 'the exon coordinate order may be wrong'
    
    # return dataframe for further analysis

    return(genome_GRCh37(df_exon_rel, root))


def genome_GRCh37(df_exon_rel, root):
    ''' Extract exome genome cordinates for build GRC37'''
# start gen pos - 5000 + data frame  nee chr start and end
    fred = root.findall('updatable_annotation/annotation_set[@type="lrg"]/mapping')
# fred = root.findall('mapping/..[@type ="lrg"]')
    print(fred)

           

    
print(xml_checker('LRG_62.xml'))
#print(xml_checker('wrong.txt'))
