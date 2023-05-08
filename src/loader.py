import time, os
from collections import defaultdict
from configs import Configs 

'''
Obtain HMMs that have sizes within the given range
'''
def obtainHMMs(indir, lower, upper, hmmbuild_paths=[]):
    Configs.log('Obtaining HMMs in given range: [{}, {}]'.format(lower, upper))
    start = time.time()

    if len(hmmbuild_paths) == 0:
        hmmbuild_paths = os.popen('find {} -name hmmbuild.model.* -type f'.format(
            indir)).read().split('\n')[:-1]

    #hmms_in_range = []
    hmm_indexes = []
    index_to_hmms = {}

    for path in hmmbuild_paths:
        lines = os.popen('head -n 15 {}'.format(path)).read().split('\n')[:-1]
        split_lines = [l.split() for l in lines]
        _dict = {x[0]: ','.join(x[1:]) for x in split_lines}
        assert 'NSEQ' in _dict, 'NSEQ not found in HMMBuild file {}!'.format(path)
        num_seq = int(_dict['NSEQ'])

        # one exception is that the backbone does not have enough sequences,
        # so only one HMM is built
        if (num_seq <= upper and num_seq >= lower) or len(hmmbuild_paths) == 1:
            #hmms_in_range.append(path)
            dirname = os.path.dirname(path)
            
            hmm_ind = int(dirname.split('/')[-1].split('_')[-1])
            hmm_indexes.append((hmm_ind, num_seq))
            index_to_hmms[hmm_ind] = (dirname, num_seq)
    
    hmm_indexes = sorted(hmm_indexes, key=lambda x: x[1], reverse=True)
    time_filter = time.time() - start
    Configs.log('Done obtaining HMMs: {}'.format(hmm_indexes))
    Configs.runtime('Time to filter HMMs based on sizes (s): {}'.format(
        time_filter))
    return hmm_indexes, index_to_hmms

'''
Obtain hmmsearch results
'''
def getHMMSearchResults(index_to_hmms):
    Configs.log('Getting all targeted (within range) HMMSearch results')
    start = time.time()
    
    scores = defaultdict(list)
    for key, val in index_to_hmms.items():
        hmm_ind = key
        hmmdir, num_seq = val

        hmmsearch_paths = os.popen('find {} -name hmmsearch.results.* -type f'.format(
            hmmdir)).read().split('\n')[:-1]

        for path in hmmsearch_paths:
            with open(path, 'r') as f:
                this_scores = eval(f.read())
                for taxon, score in this_scores.items():
                    # score[1] refers to bit-score
                    scores[taxon].append((hmm_ind, score[1]))

    # sort by bit-scores
    for taxon in scores.keys():
        scores[taxon] = sorted(scores[taxon], key=lambda x: x[1], reverse=True)

    time_get_and_rank = time.time() - start
    Configs.log('Done getting all targeted (within range) HMMSearch results')
    Configs.runtime('Time to obtain HMMSearch results and rank them (s): {}'.format(
        time_get_and_rank))
    return scores

'''
Assign query sequences to target sub-alignments based on their ranked
bit-scores. 
'''
def assignQueryToSubset(scores, hmm_indexes, index_to_hmms):
    Configs.log('Assigning query sequences to target sub-alignments')
    start = time.time()

    query_assignment = defaultdict(list)
    assigned_hmms = set()

    for taxon,score in scores.items():
        # the HMM index with the highest bit-score
        hmm_ind = score[0][0]
        query_assignment[hmm_ind].append(taxon)
        # add the hmm/sub-alignment directory to assigned_hmms
        assigned_hmms.add((hmm_ind, index_to_hmms[hmm_ind][0]))
    
    time_assign = time.time() - start
    Configs.log('Done assigning query sequences to target sub-alignments')
    Configs.runtime('Time to assign queries with the highest bit-scores (s): {}'.format(
        time_assign))
    return query_assignment, list(assigned_hmms)
