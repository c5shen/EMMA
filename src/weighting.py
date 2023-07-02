'''
Created on 10.28.2021 by Chengze Shen

Bitscore to weight calculation.
'''

import os
import time
import numpy as np
from configs import Configs

class Weights(object):
    weights = dict()
    weights_map = dict()
    ranked_bitscores = dict()
    def __init__(self):
        pass

'''
Function to read in weights from local given the taxon name
'''
def readWeights(taxon):
    infile = Configs.outdir + '/weights/w_{}.txt'.format(taxon)
    if not os.path.isfile(infile):
        return None, None
    else:
        weights, weights_map = [], [] 
        with open(infile, 'r') as f:
            line = f.read()
            taxon, raw = line.split(':')
            weights = [eval(x) for x in raw.split(';')]
            weights_map = {ind: w for (ind, w) in weights}
        return weights, weights_map

'''
Function to read in bitscores from local given the taxon name
'''
def readBitscores(taxon):
    infile = Configs.outdir + '/bitscores/b_{}.txt'.format(taxon)
    if not os.path.isfile(infile):
        return None, None
    else:
        bitscores = [] 
        with open(infile, 'r') as f:
            line = f.read()
            taxon, raw = line.split(':')
            bitscores = [eval(x) for x in raw.split(';')]
        return bitscores

'''
Function to calculate the HMM weighting, given the bitscores and sizes
of the HMMs (for a given query taxon)
inputs: ensemble of HMMs H (with their bitscores and sizes)
outputs: weights for HMMs H
'''
def calculateWeights(packed_data):
    taxon, indexes, bitscores, sizes = packed_data
    #logging.debug('working with: {}'.format(taxon))
    weights = {}
    
    assert len(indexes) == len(bitscores) == len(sizes)
    for i in range(len(bitscores)):
        score_i, size_i = bitscores[i], sizes[i]
        exponents = np.array(bitscores) - score_i \
                + np.log2(np.array(sizes) / size_i)
        denominator = np.sum(np.power(2, exponents))
        weights[indexes[i]] = 1. / denominator
    
    #num_to_retain = min(Configs.num_hmms, len(weights))
    # unlike WITCH, retain just some arbitrary number (e.g., 5) of weights
    # since we are only using THE HIGHEST ONE for assignment
    num_to_retain = min(5, len(weights))
    sorted_weights = sorted([(ind, w) for ind, w in weights.items()],
            key = lambda x: x[1], reverse=True)[:num_to_retain]
    return {taxon: tuple(sorted_weights)}

    ## write weights to local (only top k ones)
    #sorted_weights = [str(x) for x in sorted_weights]
    #with open(Configs.outdir + '/weights/w_{}.txt'.format(taxon), 'w') as f:
    #    f.write(taxon + ':' + ';'.join(sorted_weights) + '\n')
    #return None

'''
Function to write a single taxon with its ranked bitscore to local
'''
def writeQueryBitscores(packed_data):
    taxon, sorted_scores = packed_data
    str_sorted_scores = [str(x) for x in sorted_scores]

    with open(Configs.outdir + '/bitscores/b_{}.txt'.format(taxon), 'w') as f:
        f.write(taxon + ':' + ';'.join(str_sorted_scores) + '\n')
    return None

'''
Obtain and write weights to local based on bitscores
'''
def writeWeights(index_to_hmm, ranked_bitscores, pool):
    s2 = time.time()
    Configs.warning('Starting to calculate weights...')

    # - get sizes of each HMM
    all_sizes = {}
    for index, subset in index_to_hmm.items():
        # subset fields -- (dirname, alignment_ind, hmm_ind, num_seq)
        all_sizes[(subset[1], subset[2])] = subset[3]

    # iterate through each query taxon
    # write to local for each taxon and its weights
    weights, weights_map = {}, {}
    args = []
    for taxon, sorted_scores in ranked_bitscores.items():
        # each element in sorted_scores fields --
        #               (alignment_ind, hmm_ind, score[1])
        indexes = [x[0] for x in sorted_scores]
        bitscores = [x[2] for x in sorted_scores]
        sizes = [all_sizes[(x[0], x[1])] for x in sorted_scores]
        args.append((taxon, indexes, bitscores, sizes))

    all_taxon_to_weights = list(pool.map(calculateWeights, args))
    taxon_to_weights = {}
    for item in all_taxon_to_weights:
        taxon_to_weights.update(item)

    time_obtain_weights = time.time() - s2
    Configs.warning('Finished calculating weights!')
    Configs.runtime('Time to obtain weights given bitscores (s): {}'.format(
        time_obtain_weights))
    return taxon_to_weights

'''
Write weights to local as [outdir]/weights.txt
'''
def writeWeightsToLocal(taxon_to_weights, path):
    with open(path, 'w') as f:
        for taxon, weights in taxon_to_weights.items():
            f.write('{}:{}\n'.format(taxon, weights))
