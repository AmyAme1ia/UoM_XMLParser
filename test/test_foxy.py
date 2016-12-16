'''
test suite for FoxyParser.py
Uses pytest - install with: pip install -U pytest
To run tests: pytest
'''

import pandas as pd
import xml.etree.ElementTree as ET
from pandas.util.testing import assert_frame_equal
import filecmp
import os
import shutil

# set truthset variables for testing based on example LRG
FILENAME = 'LRG_files/LRG_517.xml'
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
OUTDIR = './LRG_517_output'
TESTFILE = 'LRG_517_t1.tsv'

def test_usrinput1():
    assert len(FILENAME) >1

def test_PathCheck():
    assert FILENAME.endswith('.xml')
    tree = ET.parse(FILENAME)
    root = tree.getroot()
    assert root.tag == 'lrg'
    
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

    df_gen_build.Chr = pd.to_numeric(df_gen_build.Chr)
    df_gen_build.g_start = pd.to_numeric(df_gen_build.g_start)
    df_gen_build.g_end = pd.to_numeric(df_gen_build.g_end)
    df_gen_build.strand = pd.to_numeric(df_gen_build.strand)
    assert_frame_equal(df_gen_build, DF3)
    
def test_leg():
    df_gen_build = DF3
    df = DF2
    for i in range(len(df_gen_build.Build)):
        # checks that the genome build is canonical
        if 'assembly' in str(df_gen_build.type.loc[i]):
            # check the stand orientation
            
            if str(df_gen_build.strand.loc[i]) == "-1":
                # generate a list of genomic lrg start and end position
                g_loc = df_gen_build.at[i,'g_end']
                lrg_loc_s = []
                lrg_loc_e = []
                # g_loc_e = df_gen_build.at[i,'g_start']
                
                # populate list of lrg positions
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
                
                # generate a list of lrg star positions and a ver for genomic end possition 
                g_loc = df_gen_build.at[i,'g_start']
                lrg_loc_s = []
                lrg_loc_e = []
                
                # populate list of lrg positions
                for l in range(len(df.exon_no)):
                    lrg_loc_s.append(df.start.loc[l])  
                    lrg_loc_e.append(df.end.loc[l])
                    
                # loop through calculate genomic start pos for rev strand
                exon_pos_s = [int(g_loc) + int(lrg_loc_s[x]) for x in range(len(lrg_loc_s))]
                df[(df_gen_build.Build.loc[i])+'_start'] = exon_pos_s
                                # loop through calculate genomic pos for rev strand
                exon_pos_e = [int(g_loc) + int(lrg_loc_e[x]) - 1 for x in range(len(lrg_loc_s))]
                df[(df_gen_build.Build.loc[i])+'_end'] = exon_pos_e
               
            else:
                print("Problem! DNA should only have two strands, this has more, so cant be DNA")    
    assert_frame_equal(df, DF4)