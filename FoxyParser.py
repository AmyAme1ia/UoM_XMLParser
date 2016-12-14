# Python 3.5
import xml.etree.ElementTree as ET
import pandas as pd
import warnings
warnings.filterwarnings('ignore')
import os

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
    
def loop_transcripts(root):
    ''' find the number of transcript for the record to allow iteration '''
    transcripts = root.findall('./fixed_annotation/transcript')
    print('There are',len(transcripts),'transcript(s) for this record.')
    return [transcripts[i].attrib['name'] for i in range(len(transcripts))]

def get_summary_data(root):
    ''' extract basic data from xml and store in a pandas dataframe (and gather other variables) '''    
    # find LRG id for gene
    lrg_id = root.findall('./fixed_annotation/id')[0].text
    # find gene symbol
    symbol = root.findall('./updatable_annotation/annotation_set/features/gene/symbol')[0].attrib['name']
    # find chromosome
    chromosome  = 'chr'+root.findall('./updatable_annotation/annotation_set/mapping')[0].attrib['other_name']
    # get information about strand
    strand_root = root.findall('./updatable_annotation/annotation_set/mapping/mapping_span')
    strand = strand_root[0].attrib['strand']
    return lrg_id, symbol, chromosome, strand

def get_data(root,transcript):
    ''' extract data from xml and store in a pandas dataframe (and gather other variables) '''    
    # find LRG id for gene
    lrg_id = root.findall('./fixed_annotation/id')[0].text
    # find gene symbol
    symbol = root.findall('./updatable_annotation/annotation_set/features/gene/symbol')[0].attrib['name']
    
    # define an empty dataframe to accept data
    df = pd.DataFrame(columns=['exon_no','start','end'])
    #define empty lists as an intermediary data store
    ex_no = []
    ex_start = []
    ex_end = []
    
    # parse exon data from xml file
    for item in transcript.findall('./*[@label]'):
        ex_no.append(item.attrib['label'])
        for record in item.findall('./*[1]'):
            ex_start.append(int(record.attrib['start']))
            ex_end.append(int(record.attrib['end']))
    # populate dataframe from lists
    for i in range(len(ex_no)):
        df.loc[df.shape[0]] = [ex_no[i],ex_start[i],ex_end[i]]
    return lrg_id, symbol, df

def add_sequence(df,root):
    ''' find genomic sequence and length for each exon '''
    # find genomic sequence by LRG coordinates
    genomic_sequence = root.findall('./fixed_annotation/sequence')[0].text.upper()
    # check that the genomic sequence conforms to standard DNA bases
    assert set(genomic_sequence) == set(['A', 'C', 'T', 'G']), 'Unexpected characters found in genomic sequence.'
    # add temporary indexing columns to df for sequence slice
    df['int_start'] = df.start.astype(int)
    df['int_end'] = df.end.astype(int)
    # calulate exon length and add to dataframe
    df['exon_length'] = df['end'] - df['start']
    df['seq'] = [(genomic_sequence[(df.int_start.loc[i]):(df.int_end.loc[i])]) for i in range(len(df.start))]
    # remove intermediary indexing columns
    del df['int_start']
    del df['int_end']
    # check that sequence legth matches the exon length
    for i in range(len(df.start)):
        len(df.seq.loc[i]) == df.exon_length.loc[i],
        "Sequence length doesn't match exon length"
    return df

def genome_loc(df, root):
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

    return df_gen_build

def leg (df_gen_build, df):
    '''Location of Exome in Genome'''
    
    for i in range(len(df_gen_build.Build)):
        # checks that the genome build is canonical
        if 'assembly' in str(df_gen_build.type.loc[i]):
            # check the stand orientation
            
            if str(df_gen_build.strand.loc[i]) == "-1":
                print('on reverse strand')
                # generate a list of genomic lrg start and end position
                g_loc = df_gen_build.at[i,'g_end']
                lrg_loc_s = []
                lrg_loc_e = []
                # g_loc_e = df_gen_build.at[i,'g_start']
                
                # populate list of lrg possitions
                for l in range(len(df.exon_no)):
                    lrg_loc_s.append(df.start.loc[l])    
                    lrg_loc_e.append(df.end.loc[l])
                # loop through calculate genomic start pos for rev strand
                exon_pos_s = [int(g_loc) - int(lrg_loc_s[x]) for x in range(len(lrg_loc_s))]
                df[(df_gen_build.Build.loc[i])+'_start'] = exon_pos_s
                
                # loop through calculate genomic end pos for rev strand
                exon_pos_e = [int(g_loc) - int(lrg_loc_e[x]) + 1 for x in range(len(lrg_loc_s))]
                df[(df_gen_build.Build.loc[i])+'_end'] = exon_pos_e
            
            elif str(df_gen_build.strand.loc[i]) == "1":
                print('on Forward strand')
                
                # generate a list of lrg star positions and a ver for genomic end possition 
                g_loc = df_gen_build.at[i,'g_start']
                lrg_loc_s = []
                lrg_loc_e = []
                
                # populate list of lrg possitions
                for l in range(len(df.exon_no)):
                    lrg_loc_s.append(df.start.loc[l])  
                    lrg_loc_e.append(df.end.loc[l])
                    
                # loop through calculate genomic start pos for rev strand
                exon_pos_s = [int(g_loc) + int(lrg_loc_s[x]) for x in range(len(lrg_loc_s))]
                df[(df_gen_build.Build.loc[i])+'_start'] = exon_pos_s
                                # loop through calculate genomic pos for rev strand
                exon_pos_e = [int(g_loc) + int(lrg_loc_e[x]) - 1 for x in range(len(lrg_loc_s))]
                df[(df_gen_build.Build.loc[i])+'_end'] = exon_pos_e
                
                print('genLoc:', df_gen_build.Build.loc[i])
               
            else:
                print("Problem! DNA should only have two strands, this has more, so cant be DNA")
    return df

def output_to_file(name_base,df,t,lrg_id,symbol,chromosome,strand): # from main_looper
    ''' create an output directory and writes output files to it '''
    # find path to current directory
    current_dir = os.path.dirname(os.path.realpath(__file__))
    # name new folder to contain output files based on LRG id
    new_dir_name = name_base+'_output'
    output_filename = name_base+'_'+t+'.tsv' # from main_looper
    
    # create header for output files containing basic summary data - dummy headers to allow concatenation with main dataframe
    df_head = pd.DataFrame(columns=['exon_no','start','end','exon_length','GRCh37.p13_start','GRCh37.p13_end','GRCh38.p7_start','GRCh38.p7_end','seq'])
    
    # dictionary to translate strand coding from +/-1 to +/- n header
    strand_dict = {'1':'+','-1':'-'}
    # add a series of rows containing header information and empty columns to match the size of the main dataframe
    df_head.loc[len(df_head.exon_no)] = ['#','LRG ID :',lrg_id,'','','','','','']
    df_head.loc[len(df_head.exon_no)] = ['#','Gene Symbol :',symbol,'','','','','','']
    df_head.loc[len(df_head.exon_no)] = ['#','Chromosome :',chromosome,'','','','','','']
    df_head.loc[len(df_head.exon_no)] = ['#','Strand :',strand_dict[strand],'','','','','','']
    df_head.loc[len(df_head.exon_no)] = ['#','Transcript number :',t,'','','','','','']
    df_head.loc[len(df_head.exon_no)] = ['#','','','','','','','','']
    df_head.loc[len(df_head.exon_no)] = ['# exon_no','start','end','exon_length','GRCh37.p13_start','GRCh37.p13_end','GRCh38.p7_start','GRCh38.p7_end','seq']
    
    # re-oroder main dataframe columns for output
    df = df[['exon_no','start','end','exon_length','GRCh37.p13_start','GRCh37.p13_end','GRCh38.p7_start','GRCh38.p7_end','seq']]
    # concatenate header and main dataframes prior to writing to file
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    df['exon_length'] = df['exon_length'].astype(int)
    df = pd.concat([df_head,df],axis=0)
    
    # check if new folder already exists, and create it if it doesn't
    if not os.path.exists(os.path.join(current_dir,new_dir_name)):
        os.makedirs(new_dir_name)
    df.to_csv(os.path.join(os.path.join(current_dir,new_dir_name),output_filename),sep='\t',index=False,header=False)
    return new_dir_name


def main(infile):
    '''launches main workflow for parsing LRG dtat from xml'''
    # run checks on file type
    checked = (xml_checker(infile)) 
    # checked = root
    # get LRG id, gene symbol and chromosome from file
    lrg_id, symbol, chromosome, strand = get_summary_data(checked)
    # this should be a list of markers, eg ['t1','t2',...,'tn']
    transcripts = loop_transcripts(checked)
    # loop over any available transcript records within the LRG file
    for t in transcripts:
        for transcript in checked.findall('./fixed_annotation/*[@name]'):
            if transcript.attrib['name'] == t:
                lrg_id, symbol, exon_data = get_data(checked,transcript)
                # add the exon lengths and sequences to the growing df
                exon_data_with_seq = add_sequence(exon_data,checked)

                genome_build = genome_loc(exon_data, checked)
                # genome_build is the df_gen_build dataframe
                exon_gen_pos = leg(genome_build, exon_data_with_seq)
                #print(exon_data)
                output_to_file(lrg_id,exon_gen_pos,t,lrg_id,symbol,chromosome,strand)
                print('Finished.  Please check the output folder') #.format(new_dir_name)
    #return df_genomic_coords,lrg_id,symbol,chromosome
    pass
    
    
    #print('LRG id :',lrg_id)
    #print('Gene symbol :',symbol)
    #print(exon_data_with_seq)
    #print()
    #print(genome_build)

    #return exon_data

main('LRG_517.xml') # NEED TO CHANGE SO NOT HARD CODED
    

