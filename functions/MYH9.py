#!/usr/bin/env python

from .default import *

#plot_category_classes = ['dominant', 'recessive', 'noncat']
condition_classes = ['dominant', 'recessive', 'noncat', 'eye', 'platelet', 'facial', 'other', 'noinfo']
plot_category_classes = ['dominant', 'recessive', 'noncat', 'platelet', 'facial']

def count_condition (condition_class, condition):
    condition = condition.lower()
    items = [x.strip() for x in condition.split('|')]

    def identify_condition_class (condition):
        condition = condition.lower()
        if condition.startswith('not provided') or condition.startswith('not specified'):
            return 'noinfo'
        elif condition == '\\n':
            return 'noinfo'
        elif any(x in condition for x in ['thrombocytopenia', 'platelet', 'may-hegglin', 'epistaxis']):
            return 'platelet'
        elif any(x in condition for x in ['facial', 'freckles', 'obesity']):
            return 'facial'
        elif any(x in condition for x in ['deafness', 'hearing']):
            if 'dominant' in condition:
                return 'dominant'
            elif 'recessive' in condition:
                return 'recessive'
            else:
                return 'noncat'
        else:
            print('classified as other:', condition)
            return 'other'

    return sum([condition_class == identify_condition_class(item) for item in items])

conflict_lookup = [x for x in condition_classes if x != 'noinfo' and x != 'other']

def is_conflict (variant_row):
    return sum([variant_row[x] > 0 for x in conflict_lookup]) > 1
