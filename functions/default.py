#!/usr/bin/env python

# evaluate phenotype
condition_classes = ['usher', 'dominant', 'recessive', 'noncat', 'others', 'noinfo']
def classify_condition (condition):
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

# evaluation function for variants
p_change_classes = ['stop', 'frameshift', 'missense', 'inframe', 'noncoding']
p_change_steps   = [1, 1, 2, 2, 2, 3]
def classify_p_change (p_change):
    if any(x == p_change for x in ['noncoding', 'stop', 'missense', 'frameshift']):
        return p_change
    elif p_change.startswith('inf'):
        return 'inframe'
    else:
        raise Exception('unable to classify or process:', p_change)

