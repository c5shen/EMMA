import time, os
import subprocess
from collections import defaultdict

from configs import Configs 
from helpers.alignment_tools import Alignment

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

            # read number of sequences in the corresponding subalignment
            # (not the one used in HMMs)
            subalignment_path = os.path.join('/'.join(dirname.split('/')[:-1]),
                    'subset.aln.fasta')
            assert os.path.exists(subalignment_path)
            #subalignment_size = int(os.popen('wc -l {}'.format(subalignment_path)).read().strip().split(' ')[0]) // 2
            #out = subprocess.Popen(['wc', '-l', subalignment_path],
            #                        stdout=subprocess.PIPE,
            #                        stderr=subprocess.STDOUT
            #                        ).communicate()[0]
            #subalignment_size = int(out.partition(b' ')[0]) // 2
            raw = subprocess.check_output(
                    'wc -l {}; exit 0'.format(subalignment_path),
                    stderr=subprocess.STDOUT,
                    universal_newlines=True, shell=True)
            subalignment_size = int(raw.strip().split(' ')[0]) // 2
            
            # get alignment and assignment subproblem indexes
            index = dirname.split('/')[-1]
            alignment_ind, hmm_ind = (int(x) for x in index.split('_')[1:])
            hmm_indexes.append((alignment_ind, subalignment_size, 
                hmm_ind, num_seq))
            #print(hmm_indexes[-1])
            index_to_hmms[index] = (dirname, alignment_ind, hmm_ind, num_seq)
    
    # sort by the subalignment size
    hmm_indexes = sorted(hmm_indexes, key=lambda x: x[1], reverse=True)
    time_filter = time.time() - start
    Configs.log('Done obtaining HMMs: {}'.format(hmm_indexes))
    Configs.runtime('Time to obtain HMMs based on sizes (s): {}'.format(
        time_filter))
    return hmm_indexes, index_to_hmms

'''
Obtain hmmsearch results
'''
def getHMMSearchResults(index_to_hmms):
    Configs.log('Getting all targeted (within range) HMMSearch results')
    start = time.time()
    
    scores = defaultdict(list)
    for index, val in index_to_hmms.items():
        hmmdir, alignment_ind, hmm_ind, num_seq = val

        #hmmsearch_paths = os.popen('find {} -name hmmsearch.results.* -type f'.format(
        #    hmmdir)).read().split('\n')[:-1]
        raw = subprocess.check_output(
                'find {} -name hmmsearch.results.* -type f'.format(hmmdir),
                stderr= subprocess.STDOUT,
                universal_newlines=True, shell=True)
        hmmsearch_paths = raw.strip().split('\n')

        for path in hmmsearch_paths:
            with open(path, 'r') as f:
                this_scores = eval(f.read())
                for taxon, score in this_scores.items():
                    # score[1] refers to bit-score
                    # we directly refer to the corresponding alignment subproblem
                    #scores[taxon].append((hmm_ind, score[1]))
                    scores[taxon].append((alignment_ind, hmm_ind, score[1]))

    # sort by bit-scores
    for taxon in scores.keys():
        scores[taxon] = sorted(scores[taxon], key=lambda x: x[2], reverse=True)

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

    for taxon, score in scores.items():
        # the alignment subproblem index with the highest bitscore
        # score fields -- (alignment_ind, hmm_ind, [adjusted] bitscore)
        alignment_ind = score[0][0]
        query_assignment[alignment_ind].append(taxon)
        ## add the sub-alignment directory to assigned_hmms
        #assigned_hmms.add((hmm_ind, index_to_hmms[hmm_ind][0]))
    
    time_assign = time.time() - start
    Configs.log('Done assigning query sequences to target sub-alignments')
    Configs.runtime('Time to assign queries with the highest bit-scores (s): {}'.format(
        time_assign))
    return query_assignment#, list(assigned_hmms)
'''
Function to create a local copy of the backbone alignment in upper-cases
'''
def writeTempBackbone(outdir, backbone_path):
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    tmp_backbone_path = os.path.join(outdir, 'backbone.aln.fasta')
    Configs.log('Creating a local copy of backbone alignment ' + \
            '(all letters to upper-cases) at: ' + \
            tmp_backbone_path)
    alignment = Alignment(); alignment.read_file_object(backbone_path)
    backbone_length = alignment.sequence_length()

    assert backbone_length != None, \
            'The input backbone {} is not aligned!'.format(backbone_path)

    for key in alignment.keys():
        alignment[key] = alignment[key].upper()
    alignment.write(tmp_backbone_path, 'FASTA')
    del alignment

    return tmp_backbone_path, backbone_length
