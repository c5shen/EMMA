import os, sys, time, psutil, shutil
from configs import *

from src.algorithm import DecompositionAlgorithm, SearchAlgorithm
from src.backbone import BackboneJob
from src.loader import obtainHMMs, getHMMSearchResults, \
        assignQueryToSubset
from src.writer import writeSubAlignment, writeSubQueries
from src.aligner import alignSubQueries
from src.merger import mergeAlignments

#from helpers.alignment_tools import Alignment

import multiprocessing as mp
from multiprocessing import Lock, Manager
from concurrent.futures.process import ProcessPoolExecutor
from functools import partial

# max system recursion limit hard encoding to a large number
# a temp fix for dendropy tree recursion issues
sys.setrecursionlimit(10000)

'''
Delete all unnecessary intermediate files
'''
def clearTempFiles():
    # make an empty dir for rsync removal
    blank_dir = os.path.join(Configs.outdir, 'blank')
    if not os.path.isdir(blank_dir):
        os.makedirs(blank_dir)

    dirs_to_remove = ['tree_decomp/fragment_chunks', 'tree_decomp/root']
    for _d in dirs_to_remove:
        if os.path.isdir('{}/{}'.format(Configs.outdir, _d)):
            os.system('rsync -a --delete {}/ {}/{}/'.format(blank_dir,
                Configs.outdir, _d))
            os.system('rmdir {}/{}'.format(Configs.outdir, _d))

    if os.path.isdir(blank_dir):
        os.system('rmdir {}'.format(blank_dir))


'''
initialize the pool
'''
def initiate_pool(args):
    buildConfigs(args)

def mainAlignmentProcess(args):
    m = Manager()
    lock = m.Lock()

    ######### the following codes are adapted from WITCH ##############
    # initialize the main pool at the start so that it's not memory
    # intensive
    Configs.warning('Initializing ProcessorPoolExecutor instance...')
    pool = ProcessPoolExecutor(Configs.num_cpus,
            initializer=initiate_pool, initargs=(args,))

    # 0) obtain the backbone alignment/tree and eHMMs
    # If no UPP eHMM directory provided, decompose from the backbone
    if not Configs.hmmdir:
        # if both backbone alignment/tree present
        if Configs.backbone_path and Configs.backbone_tree_path:
            pass
        else:
            # if missing backbone alignment, first perform an UPP-like
            # sequence split into backbone/query sets (i.e., randomly select
            # up to 1,000 sequences from 25% of the median length to be
            # backbone sequences)
            # Then, align the backbone sequences using MAGUS/PASTA
            # Finally, generate a FastTree2 backbone tree

            # jobs will only run if the corresponding paths are missing
            print('\nPerforming backbone alignment and/or tree estimation...')
            bb_job = BackboneJob()
            bb_job.setup()

            Configs.backbone_path, Configs.query_path = \
                    bb_job.run_alignment()
            Configs.backbone_tree_path = bb_job.run_tree()

        # after obtaining backbone alignment/tree, perform decomposition
        # and HMMSearches
        print('\nDecomposing the backbone tree...')
        decomp = DecompositionAlgorithm(Configs.backbone_path,
                Configs.backbone_tree_path, Configs.lower)
        # only run hmmbuilds on sub-alignments in ranges
        hmmbuild_paths = decomp.decomposition(Configs.lower, Configs.upper,
                lock, pool)
        print('\nPerforming targeted HMMSearches ' \
                'between the eHMM and queries...')
        search = SearchAlgorithm(hmmbuild_paths)
        hmmsearch_paths = search.search(lock, pool)
        
        # default to <outdir>/tree_comp/root
        Configs.hmmdir = Configs.outdir + '/tree_decomp/root'
    ############## codes from WITCH end ###############################

    
    # obtain HMMSearch results
    # hmm_indexes are sorted by their nums of sequences in desending order
    # index to hmms keys=(subset dirname, num seq of the sub-alignment)
    print('\nRunning the eMAFFTadd algorithm...')
    hmm_indexes, index_to_hmms = obtainHMMs(Configs.hmmdir,
            Configs.lower, Configs.upper)
    if not Configs.continue_run:
        scores = getHMMSearchResults(index_to_hmms)

        # assign queries to sub-alignments based on ranked bit-scores
        query_assignment, assigned_hmms = assignQueryToSubset(scores,
                hmm_indexes, index_to_hmms)

        # write sub-alignments to local (for running with mafft-linsi-add)
        writeSubAlignment(Configs.outdir, assigned_hmms)

        # write queries to local, split to equal-size subsets if necessary
        query_paths = writeSubQueries(Configs.outdir, hmm_indexes,
                index_to_hmms, Configs.query_path, query_assignment,
                Configs.subproblem_size)
    else:
        print('\nContinuing from the previous run...')
        query_paths = os.popen('ls {}/queries'.format(
            Configs.outdir)).read().split('\n')[:-1]
        query_paths = [Configs.outdir + '/queries/{}'.format(x) for x in
                query_paths]

    # run MAFFT-linsi-add for each subproblem
    subalignment_paths = alignSubQueries(Configs.outdir, query_paths,
            hmm_indexes)

    # merge all sub-alignments to form the final alignment
    print('\nMerging sub-alignments with transitivity...')
    mergeAlignments(Configs.output_path, subalignment_paths,
            Configs.backbone_path)

    Configs.warning('Closing ProcessPoolExecutor instance...')
    pool.shutdown()
    Configs.warning('ProcessPoolExecutor instance closed.')

    # clean temporary files (mainly decomposition files)
    clearTempFiles()

    print('\nAll done!')
