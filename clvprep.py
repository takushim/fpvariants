#!/usr/bin/env python

import sys, re, argparse
import pandas as pd
from pathlib import Path

# default values
input_filename = None
output_filename = None
output_suffix = "_prep.csv"

# parse arguments
parser = argparse.ArgumentParser(description='preprocess a clinvar variant file', \
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-o', '--output-file', default = output_filename, \
                    help='output filename ([basename]{0} if not specified)'.format(output_suffix))
parser.add_argument('input_file', default = input_filename, \
                    help='input GenBank file. Needs to include exon information.')
args = parser.parse_args()

# set arguments
def set_default_filename (filename, suffix):
    if filename is None:
        filename = str(Path(input_filename).with_suffix('')) + suffix
    return filename

input_filename = args.input_file
output_filename = set_default_filename(args.output_file, output_suffix)

# load a clinvar table
variant_table = pd.read_csv(input_filename, sep = '\t')
variant_table = variant_table[variant_table['Name'].str.match('^NM_')]

# decompose the name column
pattern = re.compile('^(NM_[\d\.]+)\((\w+)\):(c\.[\S]+)\ ?(\((p.\w+)\))?')
variant_table['refseq'] = variant_table.apply(lambda row: re.match(pattern, row['Name']).group(1), axis = 1)
variant_table['variant'] = variant_table.apply(lambda row: re.match(pattern, row['Name']).group(3), axis = 1)
variant_table['p_change'] = variant_table.apply(lambda row: re.match(pattern, row['Name']).group(5), axis = 1)

print("used refseq:", pd.factorize(variant_table['refseq'])[1].values)

# identify the position against the refseq
pattern = re.compile('^c\.(\-?[\d]+)')
def decode_variant (variant):
    match = re.match(pattern, variant)
    if match is None:
        match = re.match('^c\.(\*[\d]+)', variant)
    if match is None:
        match = re.match('^c\.[\D]+([\d]+)', variant)
    if match is None:
        raise Exception('cannot decode:', variant)
    return match.group(1)

variant_table['v_origin'] = variant_table.apply(lambda row: decode_variant(row['variant']), axis = 1)

# classify protein change
pattern = re.compile('p\.[\w]{3}[\d]+[\w]{3}')
def classify_protein_change (p_change):
    if p_change is None or p_change == '':
        return 'noncoding'
    elif 'Ter' in p_change:
        return 'stop'
    elif 'fs' in p_change:
        return 'frameshift'
    elif 'del' in p_change and 'ins' in p_change:
        return 'infrep'
    elif 'del' in p_change:
        return 'infdel'
    elif 'ins' in p_change:
        return 'infins'
    elif 'dup' in p_change:
        return 'infdup'
    elif re.match(pattern, p_change):
        return 'missense'
    else:
        raise Exception('cannot classify a protein change:', p_change)

variant_table['p_class'] = variant_table.apply(lambda row: classify_protein_change(row['p_change']), axis = 1)

# keep pathogenic or likely pathogenic variants
pattern = re.compile('^(Pathogenic|Likely pathogenic)')
variant_table = variant_table[variant_table['Clinical significance (Last reviewed)'].str.match(pattern)]
variant_table['pathogenicity'] = variant_table.apply(lambda row: re.match(pattern, row['Clinical significance (Last reviewed)']).group(1), axis = 1)

print("keeping pathogenic or likely pathogenic variants")

# output table
variant_table.to_csv(output_filename, index = False)

sys.exit()

for index, row in variant_table.iterrows():
    print(re.match(pattern, row['Name']))


