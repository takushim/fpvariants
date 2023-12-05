#!/usr/bin/env python

import sys, re, argparse
import numpy as np
import pandas as pd
from pathlib import Path
from importlib import import_module
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
eval_module_package = "functions"
eval_module_name = "default"

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
parser = argparse.ArgumentParser(description='Map variants on a coding sequence', \
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-o', '--output-file', default = output_filename, \
                    help='output filename ([basename]{0} if not specified)'.format(output_suffix))
parser.add_argument('-d', '--domain-file', default =domain_filename, \
                    help='domain filename ([basename]{0} if not specified)'.format(domain_suffix))
parser.add_argument('-v', '--variant-file', default = variant_filename, \
                    help='variant filename ([basename]{0} if not specified)'.format(variant_suffix))
parser.add_argument('-e', '--eval-module', default = eval_module_name, \
                    help='module containing evaluation methods ({0} if not specified)'.format(eval_module_name))
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
eval_module_name = args.eval_module

# Load the cDNA sequence and find the start of cds
print("loading a sequence:", input_filename)
seq_record = [record for record in SeqIO.parse(input_filename, "genbank")][0]
exons = [feature for feature in seq_record.features if feature.type == 'exon']
cdss = [feature for feature in seq_record.features if feature.type == 'CDS']
cds_offset = cdss[0].location.start
print("cds offset:", cds_offset)

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
def draw_bands (y_current, ranges, label, colors):
    for index, range in enumerate(ranges):
        axes.add_artist(Rectangle((range[0], y_current + (band_step - band_width) / 2),
                                  range[1] - range[0], band_width,
                                  color = colors[index % len(colors)]))
    axes.text(len(seq_record.seq) + 10, y_current + band_width / 2, label, va = 'center')

# start drawing from y = 0
y_current = 0

# draw exons
draw_bands(y_current, [(exon.location.start + 1, exon.location.end + 1) for exon in exons], 'exons', exon_colors)
y_current = y_current + band_step

# draw CDS(s)
draw_bands(y_current, [(cds.location.start + 1, cds.location.end + 1) for cds in cdss], 'cds', cds_colors)
y_current = y_current + band_step

# draw domain
print("loading domains:", domain_filename)
domain_table = pd.read_csv(domain_filename, header = None, delim_whitespace = True)
print(domain_table[[0, 1]])

domains = domain_table[1].to_list()
draw_bands(y_current, [(int(domain.partition('..')[0]), int(domain.partition('..')[2])) 
                       for domain in domains], 'domains', domain_colors)
y_current = y_current + band_step

# load pathogenic variants and drop duplicated variants
print("loading variants:", variant_filename)
variant_table = pd.read_csv(variant_filename)
variant_table['plot_v_origin'] = variant_table['v_origin'] + cds_offset

# classify phenotype and variants
print("evaluation module:", eval_module_name)
eval_module = import_module(eval_module_package + "." + eval_module_name)

# classify protein changes and diseases
variant_table['plot_p_change_class'] = variant_table.apply(lambda row: eval_module.classify_p_change(row['p_change_class']), axis = 1)
variant_table['plot_condition'] = variant_table.apply(lambda row: eval_module.classify_condition(row['conditions']), axis = 1)

# plot data
for c_index, condition in enumerate(eval_module.condition_classes):
    condition_table = variant_table[variant_table['plot_condition'] == condition]
    y_start = y_current + scatter_step * c_index

    axes.hlines(y_start, axes.get_xlim()[0], axes.get_xlim()[1], color = 'black', linewidth = 0.1)
    axes.text(len(seq_record.seq) + 10, y_start + scatter_step / 2, condition, va = 'center')

    for p_change_step, p_change_class in zip(eval_module.p_change_steps, eval_module.p_change_classes):
        plot_table = condition_table[condition_table['plot_p_change_class'] == p_change_class].copy()

        plot_table['plot_p_change_cum'] = plot_table.groupby('plot_v_origin').cumcount()
        plot_table['plot_y'] = y_start + p_change_step * scatter_step / (max(eval_module.p_change_steps) + 1) \
                               + scatter_delta * plot_table['plot_p_change_cum']

        axes.scatter(plot_table['plot_v_origin'], plot_table['plot_y'],
                     c = scatter_colors[eval_module.p_change_classes.index(p_change_class)],
                     marker = 'o', s = 4)

# dummy drawing to generate the legend
legend_handles = [Patch(color = scatter_colors[eval_module.p_change_classes.index(p_change)], label = p_change)
                  for p_change in eval_module.p_change_classes]
axes.legend(handles = legend_handles, ncol = len(eval_module.p_change_classes) // 2)

print("output graph:", output_filename)
figure.savefig(output_filename)

print(".")
