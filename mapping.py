#!/usr/bin/env python

import sys, re, argparse
import numpy as np
import pandas as pd
from pathlib import Path
from matplotlib import pyplot as plt
from matplotlib import colors as mcolors
from Bio import SeqIO

# default values
input_filename = None
domain_filename = None
domain_suffix = "_domains.txt"
variant_filename = None
variant_suffix = "_variants.csv"
output_filename = None
output_suffix = "_graph.png"

# key words
kwd_USH = ['Usher']
kwd_AR = ""
kwd_AD = ""
kwd_HL = ""

# parse arguments
parser = argparse.ArgumentParser(description='Map variants', \
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-o', '--output-file', default = output_filename, \
                    help='output filename ([basename]{0} if not specified)'.format(output_suffix))
parser.add_argument('-d', '--domain-file', default =domain_filename, \
                    help='domain filename ([basename]{0} if not specified)'.format(domain_suffix))
parser.add_argument('-v', '--variant-file', default = variant_filename, \
                    help='variant filename ([basename]{0} if not specified)'.format(variant_suffix))
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
domain_filename = set_default_filename(args.domain_file, domain_suffix)
variant_filename = set_default_filename(args.variant_file, variant_suffix)

# Load the cDNA sequence
seq_record = [record for record in SeqIO.parse(input_filename, "genbank")][0]

# prepare a canvas.
figure = plt.figure(figsize = (6, 4), dpi = 300)
axes = figure.add_subplot(111)
axes.spines['right'].set_visible(False)
axes.spines['left'].set_visible(False)
axes.spines['top'].set_visible(False)
axes.set_xlim(-100, len(seq_record.seq) + 100)
axes.set_xticks(np.arange(1, len(seq_record.seq) + 1, 1000))
axes.set_ylim(0, 1)
axes.get_yaxis().set_ticks([])

# draw features
exons = [feature for feature in seq_record.features if feature.type == 'exon']
exon_starts = [exon.location.start + 1 for exon in exons]
exon_ends = [exon.location.end + 1 for exon in exons]
axes.hlines([0.05] * len(exons), exon_starts, exon_ends, linewidth = 8, colors = ['black', 'gray'])

# draw CDS(s)
cdss = [feature for feature in seq_record.features if feature.type == 'CDS']
cds_starts = [cds.location.start + 1 for cds in cdss]
cds_ends = [cds.location.end + 1 for cds in cdss]
axes.hlines([0.1] * len(cdss), cds_starts, cds_ends, linewidth = 8, colors = ['blue', 'skyblue'])

# draw domain
domain_table = pd.read_csv(domain_filename, header = None, delim_whitespace = True)
domains = domain_table[1].to_list()
domain_starts = [int(domain.partition('..')[0]) for domain in domains]
domain_ends = [int(domain.partition('..')[2]) for domain in domains]
axes.hlines([0.15] * len(domains), domain_starts, domain_ends, linewidth = 8, colors = ['pink', 'limegreen'])

# load pathogenic variants and drop duplicated variants
variant_table = pd.read_csv(variant_filename)
variant_table = variant_table[variant_table['final_pathogenicity'].str.lower() == 'pathogenic']
variant_table = variant_table.drop_duplicates(subset = ['vep_hgvs_c'])

# classify variants
def classify_variant (consequence):
    consequence = consequence.lower()
    if 'stop_gained' in consequence:
        return 'stop'
    elif 'frameshift_variant' in consequence:
        return 'frameshift'
    elif 'missense_variant' in consequence:
        return 'missense'
    elif 'inflame_deletion' in consequence:
        return 'deletion'
    else:
        if 'splice' in consequence or 'intron' in consequence:
            return 'non_coding'
        else:
            raise Exception('unable to classify:', consequence)

variant_table['plot_class'] = variant_table.apply(lambda row: classify_variant(row['vep_consequence']), axis = 1)

# set plot x depending on the position on the cDNA
pattern = re.compile('^NM_[\d\.]+:c\.(-?[\d]+)\D*')
variant_table['plot_x'] = variant_table.apply(lambda row: int(re.match(pattern, row['vep_hgvs_c']).group(1)), axis = 1)

# shift the y position if multiple mutations occur in one place
pattern = re.compile('^NP_[\d\.]+:p\.([\w]{3})\d+([\w]{3})$')
variant_table['plot_missense'] = variant_table.apply(lambda row: ''.join(re.match(pattern, row['vep_hgvs_p']).group(1, 2)) \
                                                     if row['plot_class'] == 'missense' else np.nan, axis = 1)
variant_table = pd.concat([variant_table[variant_table['plot_class'] != 'missense'],
                           variant_table[variant_table['plot_class'] == 'missense'].drop_duplicates('plot_missense')]).sort_index()
variant_table['plot_y'] = variant_table.groupby('plot_x').cumcount()

print(mcolors.TABLEAU_COLORS)

for index, plot_class in enumerate(['stop', 'frameshift', 'missense', 'non-coding']):
    y_pos = 0.2 + 0.02 * index
    plot_table = variant_table[variant_table['plot_class'] == plot_class]
    axes.scatter(plot_table['plot_x'], plot_table['plot_y'] + y_pos,
                 c = list(mcolors.TABLEAU_COLORS)[index], marker = 'o', s = 4)

print(variant_table)

plt.show()