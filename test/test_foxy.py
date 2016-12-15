'''
test suite for FoxyParser.py
Uses pytest - install with: pip install -U pytest
To run tests: pytest
'''

import pandas as pd
import xml.etree.ElementTree as ET
from pandas.util.testing import assert_frame_equal

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

def test_xml_checker():
    assert FILENAME.endswith('.xml')

def test_check_file():
    tree = ET.parse(FILENAME)
    root = tree.getroot()
    assert root.tag == TAG
    
def test_check_public():
    assert ROOT.findall("./fixed_annotation/*")[4].tag == 'source'
    
def test_loop_transcripts():
    transcripts = ROOT.findall('./fixed_annotation/transcript')
    assert len(transcripts) == TRANSCRIPT_COUNT
    
def test_get_summary_data():
    lrg_id = ROOT.findall('./fixed_annotation/id')[0].text
    symbol = ROOT.findall('./updatable_annotation/annotation_set/features/gene/symbol')[0].attrib['name']
    chromosome  = 'chr'+ROOT.findall('./updatable_annotation/annotation_set/mapping')[0].attrib['other_name']
    strand_root = ROOT.findall('./updatable_annotation/annotation_set/mapping/mapping_span')
    strand = strand_root[0].attrib['strand']
    assert lrg_id == LRG_ID
    assert symbol == SYMBOL
    assert chromosome == CHROMOSOME
    assert strand == STRAND
    
def test_get_data():
    df = pd.DataFrame(columns=['exon_no','start','end'])
    ex_no = []
    ex_start = []
    ex_end = []
    for item in TRANSCRIPT[0].findall('./*[@label]'):
        ex_no.append(item.attrib['label'])
        for record in item.findall('./*[1]'):
            ex_start.append(int(record.attrib['start']))
            ex_end.append(int(record.attrib['end']))
    for i in range(len(ex_no)):
        df.loc[df.shape[0]] = [ex_no[i],ex_start[i],ex_end[i]]
    df.exon_no = pd.to_numeric(df.exon_no)    
    return df
    assert_frame_equal(df, DF1, check_dtype=False)
    
def test_add_sequence():
    genomic_sequence = ROOT.findall('./fixed_annotation/sequence')[0].text.upper()
    # check that the genomic sequence conforms to standard DNA bases
    assert set(genomic_sequence) == set(['A', 'C', 'T', 'G']), 'Unexpected characters found in genomic sequence.'
    # add temporary indexing columns to df for sequence slice
    df = DF1
    df['int_start'] = df.start.astype(int)
    df['int_end'] = df.end.astype(int)
    # calculate  exon length and add to dataframe
    df['exon_length'] = df['end'] - df['start']
    df['seq'] = [(genomic_sequence[(df.int_start.loc[i]):(df.int_end.loc[i])]) for i in range(len(df.start))]
    # remove intermediary indexing columns
    del df['int_start']
    del df['int_end']
    # check that sequence length matches the exon length
    for i in range(len(df.start)):
        len(df.seq.loc[i]) == df.exon_length.loc[i],
        "Sequence length doesn't match exon length"
    assert_frame_equal(df, DF2)

def test_genome_loc():
    GRCh_build = []
    GRCh_chr = []
    GRCh_start = []
    GRCh_end = []
    GRCh_strand = []
    GRCh_type = []
    df_gen_build = pd.DataFrame(columns=['Build','Chr', 'g_start','g_end', 'strand', 'type'])
    for item in ROOT.findall('updatable_annotation/annotation_set[@type="lrg"]/mapping'):
        GRCh_build.append(item.attrib['coord_system'])
        GRCh_chr.append(item.attrib['other_name'])
        GRCh_start.append(item.attrib['other_start'])
        GRCh_end.append(item.attrib['other_end'])
        GRCh_type.append(item.attrib['type'])
    for item in ROOT.findall('updatable_annotation/annotation_set[@type="lrg"]/mapping/mapping_span'):   
        GRCh_strand.append(item.attrib['strand'])
    for i in range(len(GRCh_build)):
        df_gen_build.loc[df_gen_build.shape[0]] = [GRCh_build[i], GRCh_chr[i], GRCh_start[i], GRCh_end[i],GRCh_strand[i], GRCh_type[i]]
    assert_frame_equal(df_gen_build, DF3, check_dtype=False)