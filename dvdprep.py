#!/usr/bin/env python

import sys, re, argparse
import pandas as pd
from pathlib import Path

# default values
input_filename = None
output_filename = None
output_suffix = "_prep.csv"

# parse arguments
parser = argparse.ArgumentParser(description='preprocess a deafness gene variant file', \
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
print("loading a variant table:", input_filename)
variant_table = pd.read_csv(input_filename)
variant_table = variant_table[~(variant_table['vep_hgvs_c'] == '\\N')]

# decompose the name column
pattern = re.compile('^(NM_[\d\.]+):(c\.[\S]+)')
variant_table['refseq'] = variant_table.apply(lambda row: re.match(pattern, row['vep_hgvs_c']).group(1), axis = 1)
variant_table['variant'] = variant_table.apply(lambda row: re.match(pattern, row['vep_hgvs_c']).group(2), axis = 1)

pattern = re.compile('^(NP_[\d\.]+):(p\.[\S]+)')
variant_table['p_change'] = variant_table.apply(lambda row: re.match(pattern, row['vep_hgvs_p']).group(2) \
                                                if row['vep_hgvs_p'] != '\\N' else '\\N', axis = 1)

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
missense_pattern = re.compile('p\.[\w]{3}[\d]+([\w]{3}|\?)')
silent_pattern = re.compile('p\.[\w]{3}[\d]+(\=)')
def classify_protein_change (p_change):
    if p_change is None or p_change == '' or p_change == '\\N':
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
    elif re.match(missense_pattern, p_change):
        return 'missense'
    elif re.match(silent_pattern, p_change):
        return 'silent'
    else:
        raise Exception('cannot classify a protein change:', p_change)

variant_table['p_change_class'] = variant_table.apply(lambda row: classify_protein_change(row['p_change']), axis = 1)

# keep pathogenic or likely pathogenic variants
pattern = re.compile('^(Pathogenic|Likely pathogenic)', re.IGNORECASE)
variant_table = variant_table[variant_table['final_pathogenicity'].str.match(pattern)]
variant_table['pathogenicity'] = variant_table['final_pathogenicity']

# conditions
variant_table['conditions'] = variant_table['final_disease']

# output table
print("output a preprocessed table:", output_filename)
variant_table.to_csv(output_filename, index = False)

print(".")


