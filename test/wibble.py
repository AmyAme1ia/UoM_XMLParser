import pandas as pd
import xml.etree.ElementTree as ET
from pandas.util.testing import assert_frame_equal
import filecmp
import os

# set truthset variables for testing based on example LRG
FILENAME = 'LRG_517.xml'
TAG = 'lrg'
TRANSCRIPT_COUNT = 1
LRG_ID = 'LRG_517'
SYMBOL = 'RB1'
CHROMOSOME = 'chr13'
STRAND = '1'
TREE = ET.parse(FILENAME)
ROOT = TREE.getroot()
TRANSCRIPT = ROOT.findall('./fixed_annotation/transcript')
DF1 = pd.read_csv('test/DF_1.tsv', sep='\t')
DF2 = pd.read_csv('test/DF_2.tsv', sep='\t')
DF3 = pd.read_csv('test/DF_3.tsv', sep='\t')
DF4 = pd.read_csv('test/DF_4.tsv', sep='\t')
T = 't1'
OUTDIRNAME = 'LRG_517_output'
OUTDIR = 'LRG_517_output'
TESTFILE = 'LRG_517_output/LRG_517_t1.tsv'


def test_output_to_file():
    name_base = LRG_ID
    df = DF4
    current_dir = os.path.dirname(os.path.realpath(__file__))
    new_dir_name = name_base+'_output'
    output_filename = name_base+'_'+T+'.tsv' # from main_looper
    df_head = pd.DataFrame(columns=['exon_no','start','end','exon_length','GRCh37.p13_start','GRCh37.p13_end','GRCh38.p7_start','GRCh38.p7_end','seq'])
    strand_dict = {'1':'+','-1':'-'}
    df_head.loc[len(df_head.exon_no)] = ['#','LRG ID :',LRG_ID,'','','','','','']
    df_head.loc[len(df_head.exon_no)] = ['#','Gene Symbol :',SYMBOL,'','','','','','']
    df_head.loc[len(df_head.exon_no)] = ['#','Chromosome :',CHROMOSOME,'','','','','','']
    df_head.loc[len(df_head.exon_no)] = ['#','Strand :',strand_dict[STRAND],'','','','','','']
    df_head.loc[len(df_head.exon_no)] = ['#','Transcript number :',T,'','','','','','']
    df_head.loc[len(df_head.exon_no)] = ['#','','','','','','','','']
    df_head.loc[len(df_head.exon_no)] = ['# exon_no','start','end','exon_length','GRCh37.p13_start','GRCh37.p13_end','GRCh38.p7_start','GRCh38.p7_end','seq']
    df = df[['exon_no','start','end','exon_length','GRCh37.p13_start','GRCh37.p13_end','GRCh38.p7_start','GRCh38.p7_end','seq']]
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    df['exon_length'] = df['exon_length'].astype(int)
    df = pd.concat([df_head,df],axis=0)
    
    if not os.path.exists(os.path.join(current_dir,new_dir_name)):
        os.makedirs(os.path.join(current_dir,new_dir_name))
    df.to_csv(os.path.join(os.path.join(current_dir,new_dir_name),output_filename),sep='\t',index=False,header=False)
    #assert name_base == OUTFILE
    assert new_dir_name == OUTDIR
    assert filecmp.cmp(os.path.join(os.path.join(current_dir,new_dir_name,output_filename)),TESTFILE)

test_output_to_file()