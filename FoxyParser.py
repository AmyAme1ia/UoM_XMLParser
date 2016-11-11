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
    return(check_public(root))

def check_public(root):
    ''' Checks that the LRG file is a public file. for pending files issues a warning regarding completeness'''
    if root.findall("./fixed_annotation/*")[4].tag == 'source':
        print('done checks')
        return(root)
    else:
        print ('Warning! this is a pending LRG file and may be subject to modification')
        return(root)

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
    print('done get data')

    return(df_exon_rel)

def add_sequence(df_exon_rel,root):
    ''' find genomic sequence for each exon '''
    # find genomic sequence by LRG coordinates
    genomic_sequence = root.findall('./fixed_annotation/sequence')[0].text.upper()
    # check that the genomic sequence conforms to standard DNA bases
    assert set(genomic_sequence) == set(['A', 'C', 'T', 'G']), 'Unexpected characters found in genomic sequence.'
    # add temporary indexing columns to df for sequence slice
    df_exon_rel['int_start'] = df_exon_rel.start.astype(int)
    df_exon_rel['int_end'] = df_exon_rel.end.astype(int)
    # calulate exon length and add to dataframe
    df_exon_rel['exon_length'] = df_exon_rel['int_end'] - df_exon_rel['int_start']
    df_exon_rel['seq'] = [(genomic_sequence[
        (df_exon_rel.int_start.loc[i]):(df_exon_rel.int_end.loc[i])]) for i in range(len(df_exon_rel.start))]
    # remove intermediary indexing columns
    df_exon_rel = df_exon_rel[['exon_no','start','end','exon_length','seq']]
    # check that sequence legth matches the exon length
    for i in range(len(df_exon_rel.start)):
        len(df_exon_rel.seq.loc[i]) == df_exon_rel.exon_length.loc[i],
        "Sequence length doesn't match exon length"
    return df_exon_rel

def genome_loc(df_exon_rel, root):
    ''' Extract exome genome cordinates for build GRC37'''
    # Generate list to extract inofrmation of genome build, chromosome, genomic start and stop possition and build assembly type
    GRCh_build = []
    GRCh_chr = []
    GRCh_start = []
    GRCh_end = []
    GRCh_strand = []
    GRCh_type = []
    
    # define an empty dataframe to accept genome build informtion from xml file
    df_gen_build = pd.DataFrame(columns=['Build','Chr', 'g_start','g_end', 'strand', 'type'])
    
    # loop through LRG file to pull out genomic information
    for item in root.findall('updatable_annotation/annotation_set[@type="lrg"]/mapping'):
        GRCh_build.append(item.attrib['coord_system'])
        GRCh_chr.append(item.attrib['other_name'])
        GRCh_start.append(item.attrib['other_start'])
        GRCh_end.append(item.attrib['other_end'])
        GRCh_type.append(item.attrib['type'])
    # pull in stand from next layer down mapping_span
    for item in root.findall('updatable_annotation/annotation_set[@type="lrg"]/mapping/mapping_span'):   
        GRCh_strand.append(item.attrib['strand'])
       
    # enter genome build data from lists into pandas dataframe
    for i in range(len(GRCh_build)):
        df_gen_build.loc[df_gen_build.shape[0]] = [GRCh_build[i], GRCh_chr[i], GRCh_start[i], GRCh_end[i],GRCh_strand[i], GRCh_type[i]]
    
    print('done genome build')

    return(df_gen_build)
    
    ####################################################################
    #find1 = df_gen_build.iat[0,2]
    #print("selecting 1 cell:", find1)
    #find2 = df_gen_build.at[0,'g_start'] 
    #print("selecting 1 cell2:", find2)
    # attempt to play with pandas and find genomic start location, not 100% accurate,may out by 1bp need to check online. Need to change variaibles. Once calculated add into df_exon_rel| new function?
    #maths = int(find1) - 5000 + int(df_exon_rel.at[0,'start'])
    #print(maths)
           
    # lrg_locus source="HGNC">FOXP3</lrg_locus
    ####################################################################

def main(infile):
    '''launches main workflow for parsing LRG dtat from xml'''
    # run checks on file type
    checked = (xml_checker(infile)) 
    # checked = root
    # retive relative exon data from LRG
    exon_data = get_data(checked)
    # exon_data is the df_exon_rel dataframe
    # calculate exon lengths and get exon sequences
    exon_data_with_seq = add_sequence(exon_data,checked)
    genome_build = genome_loc(exon_data, checked)
    # genome_build is the df_gen_build dataframe
    print(exon_data)
    print()
    print(exon_data_with_seq)
    print()
    print(genome_build)
    
    return exon_data

main('LRG_62.xml') # NEED TO CHANGE SO NOT HARD CODED
    

