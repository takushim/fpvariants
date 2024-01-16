#!/usr/bin/env python

import sys, re, argparse
import numpy as np
import pandas as pd
from pathlib import Path
from importlib import import_module
from matplotlib import pyplot as plt
from matplotlib import colors as mcolors
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from Bio import SeqIO

# default values
input_filename = None
domain_filename = None
domain_suffix = "_domains.txt"
variant_filename = None
variant_suffix = "_variants.csv"
output_filename = None
output_suffix = "_graph.png"
sequence_range = [0, 0]
eval_module_package = "functions"
eval_module_name = "default"

# graph parameters
band_step = 0.04
exon_width = 0.025
domain_width = 0.04
cds_width = 0.04
#cds_width = 0.005
exon_colors = ['black', 'gray']
cds_colors = ['gray', 'blue'] # the second one will not be used.
domain_colors = ['pink', 'lightgreen'] # if not specified
scatter_step = 0.1
scatter_delta = 0.01
scatter_colors = list(mcolors.TABLEAU_COLORS)
scatter_marker_default = 'o'
scatter_marker_conflict = 'X'
invert_yaxis = False

# canvas settings
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['figure.figsize'] = (4, 3)
plt.rcParams['figure.dpi'] = 300
plt.rcParams['figure.frameon'] = False
plt.rcParams['font.size'] = 8
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['axes.linewidth'] = 1
plt.rcParams['axes.spines.left'] = False
plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.spines.top'] = invert_yaxis
plt.rcParams['axes.spines.bottom'] = not invert_yaxis
plt.rcParams['xtick.top'] = invert_yaxis
plt.rcParams['xtick.bottom'] = not invert_yaxis
plt.rcParams['xtick.labeltop'] = invert_yaxis
plt.rcParams['xtick.labelbottom'] = not invert_yaxis
plt.rcParams['ytick.left'] = False
plt.rcParams['ytick.right'] = False
plt.rcParams['ytick.labelleft'] = False
plt.rcParams['ytick.labelright'] = False
plt.rcParams['lines.linewidth'] = 0.5
plt.rcParams['lines.markersize'] = 3
plt.rcParams['scatter.edgecolors'] = 'none'
plt.rcParams['legend.frameon'] = False
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['legend.loc'] = 'lower right' if invert_yaxis else 'upper right'
plt.rcParams['savefig.transparent'] = True

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
parser.add_argument('-s', '--sequence-range', default = sequence_range, nargs = 2, type = int, \
                    help='range of sequence to draw (MIN, MAX). Inf for MAX = 0.')
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
sequence_range = args.sequence_range

# Load the cDNA sequence and find the start of cds
print("loading a sequence:", input_filename)
seq_record = [record for record in SeqIO.parse(input_filename, "genbank")][0]
exons = [feature for feature in seq_record.features if feature.type == 'exon']
cdss = [feature for feature in seq_record.features if feature.type == 'CDS']
cds_offset = cdss[0].location.start
print("cds offset:", cds_offset)

# prepare a canvas.
figure = plt.figure()
axes = figure.add_subplot(111)
xlim_min = sequence_range[0] - 100
xlim_max = len(seq_record.seq) + 100 if sequence_range[1] == 0 else sequence_range[1] + 100
print("using xlim:", xlim_min, xlim_max)
axes.set_xlim(xlim_min, xlim_max)
axes.set_ylim(0, 1)
if invert_yaxis:
    axes.invert_yaxis()

# load pathogenic variants and drop duplicated variants
print("loading variants:", variant_filename)
variant_table = pd.read_csv(variant_filename)
variant_table['plot_v_origin'] = variant_table['v_origin'] + cds_offset

# classify phenotype and variants
print("evaluation module:", eval_module_name)
eval_module = import_module(eval_module_package + "." + eval_module_name)

# classify protein changes and diseases
variant_table['plot_p_change_class'] = variant_table.apply(lambda row: eval_module.identify_plot_p_change_class(row['p_change_class']), axis = 1)
variant_table['plot_marker_color'] = variant_table.apply(lambda row: scatter_colors[eval_module.plot_p_change_classes.index(row['plot_p_change_class'])], axis = 1)

# start drawing from y = 0
y_current = 0
def draw_bands (y_current, ranges, width, colors):
    for index, range in enumerate(ranges):
        axes.add_artist(Rectangle((range[0], y_current + (band_step - width) / 2), range[1] - range[0], width,
                                  color = colors[index % len(colors)], linestyle = None, linewidth = 0))

# draw exons
draw_bands(y_current, [(exon.location.start + 1, exon.location.end + 1) for exon in exons], exon_width, exon_colors)
y_current = y_current + band_step

# plot data
for c_index, plot_category in enumerate(eval_module.plot_category_classes):
    plot_category_table = variant_table[variant_table.apply(lambda row: eval_module.is_condition_to_plot(plot_category, row), axis = 1)]
    y_start = y_current + scatter_step * c_index

    for plot_step in range(max(eval_module.plot_p_change_steps) + 1):
        plot_p_change_classes = [c for c, s in zip(eval_module.plot_p_change_classes, eval_module.plot_p_change_steps) if s == plot_step]
        plot_p_change_table = plot_category_table[plot_category_table['plot_p_change_class'].isin(plot_p_change_classes)].copy()
        plot_p_change_table['plot_p_change_cum'] = plot_p_change_table.groupby('plot_v_origin').cumcount()
        plot_p_change_table['plot_y'] = y_start + plot_step * scatter_step / (max(eval_module.plot_p_change_steps) + 1) \
                                      + scatter_delta * plot_p_change_table['plot_p_change_cum']

        conflict_mask = plot_p_change_table['conflict']
        axes.scatter(plot_p_change_table['plot_v_origin'][~conflict_mask], plot_p_change_table['plot_y'][~conflict_mask], 
                     c = plot_p_change_table['plot_marker_color'][~conflict_mask], marker = scatter_marker_default)
        axes.scatter(plot_p_change_table['plot_v_origin'][conflict_mask], plot_p_change_table['plot_y'][conflict_mask], 
                     c = plot_p_change_table['plot_marker_color'][conflict_mask], marker = scatter_marker_conflict)

    axes.text(xlim_max, y_start + scatter_step / 2, "{0}: {1}".format(len(plot_category_table), plot_category), va = 'center')

# lines between plot_category_classes. omit the last line.
for c_index in range(len(eval_module.plot_category_classes) - 1):
    y_start = y_current + scatter_step * c_index
    axes.hlines(y_start + scatter_step, axes.get_xlim()[0], axes.get_xlim()[1], color = 'black', linestyles = 'dashed')

# end drawing variants
y_current = y_current + scatter_step * len(eval_module.plot_category_classes)

# draw CDS(s)
draw_bands(y_current, [(cds.location.start + 1, cds.location.end + 1) for cds in cdss], cds_width, cds_colors)

# draw domain
print("loading domains:", domain_filename)
domain_table = pd.read_csv(domain_filename, header = None, sep = '\t')
domain_table = domain_table.apply(lambda x: x.str.strip() if x.dtype == 'object' else x)
print(domain_table[[0, 1]])

def domain_range (item):
    return (int(item.partition('..')[0].replace(",","")), int(item.partition('..')[2].replace(",","")))

if len(domain_table.columns) > 5:
    domain_colors = domain_table[5].to_list()
draw_bands(y_current, [domain_range(x) for x in domain_table[1].to_list()], domain_width, domain_colors)
y_current = y_current + band_step

# legend
legend_handles = [Line2D([0], [0], color = scatter_colors[eval_module.plot_p_change_classes.index(p_change)],
                         label = p_change, marker = 'o', linestyle = 'none', linewidth = 0)
                  for p_change in eval_module.plot_p_change_classes]
axes.legend(handles = legend_handles, ncol = len(eval_module.plot_p_change_classes) // 2)

# output graph
print("output graph:", output_filename)
for o in figure.findobj():
    o.set_clip_on(False)

figure.savefig(output_filename)

print(".")
