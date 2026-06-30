import pandas as pd
import numpy as np
import re
from Bio import SeqIO
import colorcet
from matplotlib import colors
import string

# Paths

base_dir = '/home/mathieu/postdoc_heasley/long_read_project/'
ref_dir = '/home/mathieu/postdoc_heasley/long_read_project/data/ref/'
fig_path = '/home/mathieu/postdoc_heasley/long_read_project/fig/'

# Input functions

def parse_paf(path):
    with open(path) as fi:
        paf = pd.DataFrame([line.split('\t') for line in fi.read().splitlines()]).iloc[:, :17]\
        .astype({0:str,1:int,2:int,3:int,4:str,5:str,6:int,7:int,8:int,9:int,10:int,11:int})\
        .rename({0:'rid',1:'qlen',2:'qstart',3:'qend',4:'strand',
                 5:'subject',6:'slen',7:'sstart',8:'send',
                 9:'matches', 10:'aln_block_len',11:'mapq'}, axis=1)
    return paf

def paf2bed(paf_path, subject=None):
    paf = parse_paf(paf_path)
    if subject != None:
        paf = paf.loc[paf['subject'] == subject]
    paf['name'] = paf.apply(lambda x: f'{x["subject"]}_{x["sstart"]}-{x["send"]}({x["strand"]})', axis=1)
    bed = paf[['rid', 'qstart', 'qend', 'name','mapq','strand']]
    return bed

def parse_psl(path):
    with open(path) as fi:
        psl = pd.DataFrame([line.split('\t') for line in fi.read().splitlines()]).iloc[:, :21]\
        .rename({0:'matches', 1:'misMatches', 2:'repMatches', 3:'nCount',
                 4:'qNumInsert', 5:'qBaseInsert', 6:'tNumInsert', 7:'tBaseInsert', 8:'strand',
                 9:'qName', 10:'qSize', 11:'qStart', 12:'qEnd',
                 13:'tName', 14:'tSize', 15:'tStart', 16:'tEnd',
                 17:'blockCount', 18:'blockSizes', 19:'qStarts', 20:'tStarts'}, axis=1)
        psl = psl.astype(dict([(c,str) if c in ['strand', 'qName', 'tName', 'blockSizes', 'qStarts', 'tStarts'] else (c,np.int32) for c in psl.columns]))
        for c in ['blockSizes', 'qStarts', 'tStarts']:
            psl[c] = psl[c].apply(lambda x: np.array(x.strip(',').split(',')).astype(np.int32))
    return psl

def parse_gff(path, exploded=True):
    gff = pd.read_csv(path, sep='\t', comment='#', header=None)
    gff[0] = gff[0].astype(str)
    gff[1] = gff[1].astype(str)
    gff[2] = gff[2].astype(str)
    gff[3] = gff[3].astype(np.int64)
    gff[4] = gff[4].astype(np.int64)
    gff[5] = gff[5].replace({'.':np.nan}).astype(np.float64)
    gff[6] = gff[6].astype(str)
    gff[7] = gff[7].replace({'.':np.nan}).astype(np.float64)
    gff.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    
    if exploded:
        gff = pd.concat([gff.iloc[:, :8], gff['attributes'].apply(lambda x: pd.Series(dict([i.split('=') for i in x.split(';')])))], axis=1)

    return gff

# Reference genome

with open(f'{ref_dir}S288C.fasta') as fi:
    ref_genome = [seq for seq in SeqIO.parse(fi, 'fasta')]

nuclear_chromosomes = [seq.id for seq in ref_genome][:-1]
nuc_chr_len = {seq.id:len(seq.seq) for seq in ref_genome if seq.id in nuclear_chromosomes}
nuc_chr_alias = dict(zip(nuclear_chromosomes, ['chrI','chrII','chrIII','chrIV',
                                               'chrV','chrVI','chrVII','chrVIII',
                                               'chrIX','chrX','chrXI','chrXII',
                                               'chrXIII','chrXIV','chrXV','chrXVI']))
nuc_chr_alias_num = dict(zip(nuclear_chromosomes, ['chr1','chr2','chr3','chr4',
                                               'chr5','chr6','chr7','chr8',
                                               'chr9','chr10','chr11','chr12',
                                               'chr13','chr14','chr15','chr16']))

gff = pd.read_csv(f'{ref_dir}S288C.gff', sep='\t', comment='#', header=None)
centromeres = gff.loc[gff[2]=='centromere'].groupby(0).apply(lambda x: x[[3,4]].values.flatten().mean(), include_groups=False)

chr_list = pd.read_csv(f'{ref_dir}S288C.chr_list.tsv', sep='\t')
chr_ctg_to_alias = chr_list.set_index('contig_name')['alias'].to_dict()
chr_alias_to_ctg = {j:i for (i,j) in chr_ctg_to_alias.items()}

tig_off = pd.Series(dict([(seq.id, len(seq.seq)) for seq in ref_genome]), name='len')
tig_off = pd.concat([tig_off, tig_off.cumsum().rename('cumul_end')], axis=1)
tig_off['cumul_start'] = np.concatenate([[0], tig_off['len'].iloc[:-1].cumsum().values])

# Strain metadata

batches_info = pd.read_csv(f'{base_dir}script/batches_info.csv')

strains_info = pd.DataFrame(batches_info['strain'].unique(), columns=['strain'])
strains_info.index = strains_info['strain']
strains_info = strains_info.drop('JAY270')

collection = pd.read_excel('/home/mathieu/postdoc_heasley/collection/mccusker_collection_wild_annot.xlsx').set_index('ID')

strains_info = pd.merge(strains_info, collection[['clade', 'source', 'source_simplified', 'Y\' element']], left_index=True, right_index=True) 
strain_color = {s:colorcet.glasbey_category10[i] for i,s in enumerate(list(strains_info.index))}

clades_unique = pd.read_csv('/home/mathieu/postdoc_heasley/short_read_project/script/clades_unique.csv', index_col=0)
clade_order = clades_unique.index
clade_color = clades_unique['mfc'].to_dict()
clade_idx = dict(zip(clade_order, range(len(clade_order))))
strains_info['clade_idx'] = strains_info['clade'].apply(lambda x: clade_idx.get(x, np.nan))

with open(f'{base_dir}assemblies_fisher/accessions.txt') as handle:
    ena_accessions = handle.read().splitlines()

odonnell_S1 = pd.read_excel('/home/mathieu/Dropbox/Travail/postdoc_heasley/papiers/odonnell_natgen2023_SUPP.xlsx', 
                            sheet_name='Supp Table 1', skiprows=2)

ena_accession_match = re.compile(r'(CAS\w{3})01;')
odonnell_S1['ENA'] = odonnell_S1['assembly_accession'].apply(lambda x: re.match(ena_accession_match, x))
odonnell_S1['ENA'] = odonnell_S1['ENA'].apply(lambda x: x.group(1) if x != None else np.nan)

std_name_match = re.compile(r'([ABC][A-Z]{2})')
odonnell_S1['std_name'] = odonnell_S1['standardized_name'].apply(lambda x: re.match(std_name_match, x))
odonnell_S1['std_name'] = odonnell_S1['std_name'].apply(lambda x: x.group(1) if x != None else np.nan)

odonnell_S1.index = odonnell_S1['std_name'].values

ena_to_std_dict = dict(zip(odonnell_S1['ENA'], odonnell_S1['std_name']))

# Import 1011 strains metadata
strains_1011 = pd.read_csv(f'/home/mathieu/Dropbox/Travail/postdoc_heasley/short_read_project/script/1011_genomes/peter_nature2018_S1.tsv', sep='\t')\
.set_index('Standardized name', drop=False)
# Annotate clades
strains_1011['Clades_original'] = strains_1011['Clades'].fillna('Unclustered').apply(lambda x: x.strip(' '))
strains_1011['Clades'] = strains_1011['Clades'].fillna('Unclustered').apply(lambda x: x.strip(' '))
# Manually erase the subclades of Wine/European
strains_1011['Clades'] = strains_1011['Clades'].replace({f'1. Wine/European (subclade {i})':'1. Wine/European' for i in range(1,5)})
pattern_clade_no = re.compile(r'(\d+)\. \w+')
clade_no_unmatched = {'M1. Mosaic region 1': 27, 'M2. Mosaic region 2':28, 'M3. Mosaic region 3':29, 'Unclustered': 30}
strains_1011['clade_no'] = strains_1011['Clades'].apply(lambda x: int(re.match(pattern_clade_no, x).group(1)) if re.match(pattern_clade_no, x) else clade_no_unmatched[x])
strains_1011['Total number of SNPs'] = strains_1011['Total number of SNPs'].apply(lambda x: x.replace(',', '')).astype(np.int64)

clades_unique = strains_1011.groupby('clade_no').apply(lambda x: x['Clades'].iloc[0], include_groups=False).rename('Clades')
clades_unique.loc[1] = '1. Wine/European'
clades_unique.loc[30] = 'Soil/tree'
clades_unique.loc[31] = 'Californian sour figs'
clades_unique.loc[32] = 'Unclustered'

clades_unique = clades_unique.reset_index()
clades_unique.index = clades_unique['Clades'].values
# Manually populate with the new YJM strains
clades_unique['mfc'] = [c+[0.8] for c in colorcet.glasbey_hv[:clades_unique.shape[0]]]
#clades_unique.to_csv(f'{base_dir}script/clades_unique.csv')

# Plotting variables

source_simplified_dict = {'clinical':0,
                          'wine':1,
                          'fruit':2,
                          'fermented_foods':1,
                          'brewing':1,
                          'unknown':3,
                          'fermented_beverages':1,
                          'soil':2,
                          'industry':1,
                          'tree':2,
                          'sake':1,
                          'insect':2,
                          'distillery':1,
                          'dairy':1,
                          'baking':1,
                          'bark':2,
                          'human':1,
                          'palm_nectar':2,
                          'commercial':1,
                          'sewage':1, 
                          'molasses':1,
                          'natural_fermentation':1,
                          'grain':1,
                          'water':1, 
                          'evolution_canyon':2, 
                          'plants':2,
                          'rotten_wood':2,
                          'mushrooms':2,
                          'spoiled_foods':1,
                          'fermentation':1,
                          'animal':1,
                          'alpechin':1,
                          'vegetables':1,
                          'slime_flux':2,
                          'flower':2,
                          'Domesticated':1,
                          'Human':0,
                          'Wild':2,
                          'Laboratory':1}

source_simplified_order = ['clinical', 'domesticated', 'environmental', 'unknown']
source_simplified_int = dict(zip(source_simplified_order, range(4)))
source_simplified_color = dict(zip(source_simplified_order, ['red', 'royalblue', 'limegreen', '0.7']))
cmap_source_simplified_color = colors.LinearSegmentedColormap.from_list('source_simplified_color', ['red', 'royalblue', 'limegreen', '0.7'], N=4)

repeats_sequences = {}
repeats_lengths = {}

ty5_types = []
x_types = []
#yprime_types = [f'yprime.{i}.trim_cons' for i in [1,2,3,4,5,6,7,8,9,10,11,12,13,18,19,20,21]]
yprime_types = [f'yprime.{i}.trim_cons' for i in [9,10,2,5,12,4,6,20,1,8,11,3,18,7,19,13,21]]

yprime_type_pattern = re.compile(r'yprime.(\d+).trim_cons')

yprime_types_encoding = dict(zip(yprime_types, list(string.ascii_lowercase)))

for seq in SeqIO.parse(f'{ref_dir}repeats.fasta', 'fasta'):
    name = seq.id
    if name[:3] == 'TY5':
        ty5_types.append(name)
    elif name[:6] == 'X_core':
        x_types.append(name)

    repeats_sequences[name] = seq
    repeats_lengths[name] = len(seq.seq)
repeats_lengths['TG_repeat'] = 0

at_subtelomeres = ['TG_repeat'] + ty5_types + x_types + yprime_types
at_subtelomeres_woTY5 = ['TG_repeat'] + x_types + yprime_types
at_color = dict(zip(at_subtelomeres, ['dodgerblue'] + ['limegreen']*2 + ['k' for i in range(16)] + colorcet.glasbey_warm[1:]))
at_color['Y_prime_element.cons_trim'] = 'gold'
at_alias = dict(zip(at_subtelomeres, ['TG repeat', 'Ty5', 'Ty5'] + ['X' for x in x_types] +\
                    [f'Y\' type {int(re.match(yprime_type_pattern, y).group(1))}' for y in yprime_types]))

atcg_num = {'-':-1, 'n':0, 'a':1, 't':2, 'c':3, 'g':4}
atcg_cmap = colors.LinearSegmentedColormap.from_list('ATCG', ['white', '0.4', 'limegreen', 'red', 'royalblue', 'black'], N=6)