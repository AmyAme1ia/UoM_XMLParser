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
DF1 = pd.read_csv('DF_1.tsv', sep='\t')
DF2 = pd.read_csv('DF_2.tsv', sep='\t')
DF3 = pd.read_csv('DF_3.tsv', sep='\t')
DF4 = pd.read_csv('DF_4.tsv', sep='\t')



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
        
    return df
df = test_get_data()

df.exon_no = pd.to_numeric(df.exon_no)
#DF1.exon_no = pd.to_numeric(DF1.exon_no)

print('DF1')    
print(DF1.head())
print('df')
print(df.head())

print(df == DF1)
print(df.columns.to_series().groupby(df.dtypes).groups)
print(DF1.columns.to_series().groupby(df.dtypes).groups)