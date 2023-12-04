#!/usr/bin/env python

import sys, re, argparse
import numpy as np
import pandas as pd
from pathlib import Path
from matplotlib import pyplot as plt
from matplotlib import colors as mcolors
from matplotlib.patches import Rectangle, Patch
from Bio import SeqIO

# default values
input_filename = None
domain_filename = None
domain_suffix = "_domains.txt"
variant_filename = None
variant_suffix = "_variants.csv"
output_filename = None
output_suffix = "_graph.png"
mapping_offset = 0

# graph parameters
band_step = 0.05
band_width = 0.04
exon_colors = ['black', 'gray']
cds_colors = ['lightskyblue', 'green'] # the second one will not be used.
domain_colors = ['pink', 'lightgreen']

scatter_step = 0.1
scatter_delta = 0.005
scatter_colors = list(mcolors.TABLEAU_COLORS)

# parse arguments
parser = argparse.ArgumentParser(description='Map variants', \
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-o', '--output-file', default = output_filename, \
                    help='output filename ([basename]{0} if not specified)'.format(output_suffix))
parser.add_argument('-d', '--domain-file', default =domain_filename, \
                    help='domain filename ([basename]{0} if not specified)'.format(domain_suffix))
parser.add_argument('-v', '--variant-file', default = variant_filename, \
                    help='variant filename ([basename]{0} if not specified)'.format(variant_suffix))
parser.add_argument('-m', '--mapping-offset', default = mapping_offset, type = int, \
                    help='offset of variant mapping against the GenBank sequence')
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

mapping_offset = args.mapping_offset

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

# drawing function
def draw_bands (y_current, ranges, colors):
    for index, range in enumerate(ranges):
        axes.add_artist(Rectangle((range[0], y_current + (band_step - band_width) / 2),
                                  range[1] - range[0], band_width,
                                  color = colors[index % len(colors)]))

# start drawing from y = 0
y_current = 0

# draw features
exons = [feature for feature in seq_record.features if feature.type == 'exon']
draw_bands(y_current, [(exon.location.start + 1, exon.location.end + 1) for exon in exons], exon_colors)
axes.text(len(seq_record.seq) + 10, y_current + band_width / 2, "exons", va = 'center')
y_current = y_current + band_step

# draw CDS(s)
cdss = [feature for feature in seq_record.features if feature.type == 'CDS']
draw_bands(y_current, [(cds.location.start + 1, cds.location.end + 1) for cds in cdss], cds_colors)
axes.text(len(seq_record.seq) + 10, y_current + band_width / 2, "cds", va = 'center')
y_current = y_current + band_step

# draw domain
domains = pd.read_csv(domain_filename, header = None, delim_whitespace = True)[1].to_list()
draw_bands(y_current, [(int(domain.partition('..')[0]), int(domain.partition('..')[2])) 
                       for domain in domains], domain_colors)
axes.text(len(seq_record.seq) + 10, y_current + band_width / 2, "domains", va = 'center')
y_current = y_current + band_step

# load pathogenic variants and drop duplicated variants
variant_table = pd.read_csv(variant_filename)
variant_table = variant_table[variant_table['final_pathogenicity'].str.lower() == 'pathogenic']
variant_table = variant_table.drop_duplicates(subset = ['vep_hgvs_c'])

# evaluation function for phenotype
phenotype_classes = ['usher', 'hl_dom', 'hl_rec', 'hl_uncat', 'others']
def classify_phenotype (phenotype):
    phenotype = phenotype.lower()
    if 'usher' in phenotype:
        return 'usher'
    elif any(x in phenotype for x in ['deafness', 'hearing']):
        if 'retinitis' in phenotype:
            return 'usher'
        elif 'dominant' in phenotype:
            return 'hl_dom'
        elif 'recessive' in phenotype:
            return 'hl_rec'
        else:
            return 'hl_uncat'
    else:
        return 'others'

# evaluation function for variants
variant_classes = ['stop', 'frameshift', 'missense', 'deletion', 'insertion', 'noncoding']
variant_plots   = [1, 1, 2, 2, 2, 3]
def classify_variant (consequence):
    consequence = consequence.lower()
    if 'stop_gained' in consequence:
        return 'stop'
    elif 'frameshift_variant' in consequence:
        return 'frameshift'
    elif 'missense_variant' in consequence:
        return 'missense'
    elif 'protein_altering_variant' in consequence:
        return 'missense'
    elif 'inframe_deletion' in consequence:
        return 'deletion'
    elif 'inframe_insertion' in consequence:
        return 'insertion'
    else:
        if 'splice' in consequence or 'intron' in consequence:
            return 'noncoding'
        else:
            raise Exception('unable to classify:', consequence)

# classify variants
variant_table['phenotype_class'] = variant_table.apply(lambda row: classify_phenotype(row['final_disease']), axis = 1)
variant_table['variant_class'] = variant_table.apply(lambda row: classify_variant(row['vep_consequence']), axis = 1)

# check duplicates
pattern = re.compile('^NP_[\d\.]+:p\.([\w]{3})\d+([\w]{3})')
variant_table['dup_check'] = variant_table.apply(lambda row: row['phenotype_class'] + '-' + '-'.join(re.match(pattern, row['vep_hgvs_p']).group(1, 2)) \
                                                 if row['variant_class'] == 'missense' else np.nan, axis = 1)
variant_table = pd.concat([variant_table[variant_table['variant_class'] != 'missense'],
                           variant_table[variant_table['variant_class'] == 'missense'].drop_duplicates('dup_check')]).sort_index()

# set plot x depending on the position on the cDNA
pattern = re.compile('^NM_[\d\.]+:c\.(-?[\d]+)\D*')
variant_table['plot_x'] = variant_table.apply(lambda row: int(re.match(pattern, row['vep_hgvs_c']).group(1)), axis = 1)
if mapping_offset != 0:
    variant_table['plot_x'] += mapping_offset
    print("Mapping of variants are offset by {0} to the 3'-end against the GenBank sequence".format(mapping_offset))

# plot data
for p_index, phenotype in enumerate(phenotype_classes):
    phenotype_table = variant_table[variant_table['phenotype_class'] == phenotype]
    y_start = y_current + scatter_step * p_index

    axes.hlines(y_start, axes.get_xlim()[0], axes.get_xlim()[1], color = 'black', linewidth = 0.1)
    axes.text(len(seq_record.seq) + 10, y_start + scatter_step / 2, phenotype, va = 'center')

    for variant_plot, variant_class in zip(variant_plots, variant_classes):
        plot_table = phenotype_table[phenotype_table['variant_class'] == variant_class].copy()

        plot_table['plot_cum'] = plot_table.groupby('plot_x').cumcount()
        plot_table['plot_y'] = y_start + variant_plot * scatter_step / (max(variant_plots) + 1) \
                               + scatter_delta * plot_table['plot_cum']

        axes.scatter(plot_table['plot_x'], plot_table['plot_y'],
                     c = scatter_colors[variant_classes.index(variant_class)],
                     marker = 'o', s = 4)

# dummy drawing to generate the legend
legend_handles = [Patch(color = scatter_colors[variant_classes.index(variant)], label = variant)
                  for variant in variant_classes]
axes.legend(handles = legend_handles, ncol = len(variant_classes) // 2)

plt.show()
