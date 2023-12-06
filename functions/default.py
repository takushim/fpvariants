#!/usr/bin/env python

import pandas as pd

### members for pre-processing
# evaluate phenotype
condition_classes = ['usher', 'dominant', 'recessive', 'noncat', 'systemic', 'other', 'noinfo']

def count_condition_clinvar (condition_class, condition):
    condition = condition.lower()
    items = [x.strip() for x in condition.split('|')]
    return sum([condition_class == identify_condition_class(item) for item in items])

# usually the same classification works for the dvd database
count_condition_dvd = count_condition_clinvar

def identify_condition_class (condition):
    condition = condition.lower()
    if condition.startswith('not provided'):
        return 'noinfo'
    elif 'usher' in condition:
        return 'usher'
    elif any(x in condition for x in ['deafness', 'hearing']):
        if 'retinitis' in condition:
            return 'usher'
        elif 'dominant' in condition:
            return 'dominant'
        elif 'recessive' in condition:
            return 'recessive'
        else:
            return 'noncat'
    else:
        return 'others'

conflict_lookup = [x for x in condition_classes if x != 'noinfo']

def is_conflict (variant_row):
    return sum([variant_row[x] > 0 for x in conflict_lookup]) > 1

### members for plotting
# conditions to make a scatter plot
plot_category_classes = ['usher', 'nonsynd', 'dominant', 'recessive', 'noncat', 'systemic', 'other', 'noinfo']

# categorize conditions allowing duplicates
def is_condition_to_plot (plot_category, variant_row):
    if plot_category == 'nonsynd':
        if any(variant_row[x] > 0 for x in ['dominant', 'recessive', 'noncat']):
            return True
    elif variant_row[plot_category] > 0:
        return True
    return False

# protein changes to draw plots
plot_p_change_classes = ['stop', 'frameshift', 'missense', 'inframe', 'silent', 'noncoding']

# row of drawing (from bottom to top, None = skip)
plot_p_change_steps   = [1, 1, 2, 2, 2, 3]

# evaluation function for variants
def identify_plot_p_change_class (p_change):
    if any(x == p_change for x in ['noncoding', 'stop', 'missense', 'silent', 'frameshift']):
        return p_change
    elif p_change.startswith('inf'):
        return 'inframe'
    else:
        raise Exception('unable to classify:', p_change)
