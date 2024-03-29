#!/usr/bin/env python

from .default import *

plot_category_classes = ['usher', 'dominant', 'recessive', 'noncat', 'eye']

condition_classes = ['usher', 'dominant', 'recessive', 'noncat', 'eye', 'systemic', 'other', 'noinfo']

def count_condition (condition_class, condition):
    condition = condition.lower()
    items = [x.strip() for x in condition.split('|')]

    def identify_condition_class (condition):
        condition = condition.lower()
        if condition.startswith('not provided') or condition.startswith('not specified'):
            return 'noinfo'
        elif condition == '\\n':
            return 'noinfo'
        elif 'usher' in condition:
            return 'usher'
        elif any(x in condition for x in ['deafness', 'hearing']):
            if 'retinitis pigmentosa' in condition:
                return 'usher'
            elif 'retinal dystrophy' in condition:
                return 'usher'
            elif 'dominant' in condition:
                return 'dominant'
            elif 'recessive' in condition:
                return 'recessive'
            else:
                return 'noncat'
        elif any(x in condition for x in ['retinal', 'retinitis']):
            return 'eye'
        else:
            print('classified as other:', condition)
            return 'other'

    return sum([condition_class == identify_condition_class(item) for item in items])


