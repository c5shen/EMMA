import os, sys, time, psutil, shutil
from configs import *

from src.algorithm import DecompositionAlgorithm, SearchAlgorithm
from src.backbone import BackboneJob
from src.loader import obtainHMMs, getHMMSearchResults, \
        assignQueryToSubset, writeTempBackbone
from src.writer import writeSubAlignment, writeSubQueries
from src.aligner import alignSubQueries
from src.merger import mergeAlignments
from src.weighting import writeWeights, writeWeightsToLocal, \
        readWeightsFromLocal

#from helpers.alignment_tools import Alignment

import multiprocessing as mp
from multiprocessing import Lock, Manager
from concurrent.futures.process import ProcessPoolExecutor
from functools import partial

# max system recursion limit hard encoding to a large number
# a temp fix for dendropy tree recursion issues
sys.setrecursionlimit(25000)

'''
Delete all unnecessary intermediate files
'''
def clearTempFiles():
    # make an empty dir for rsync removal
    blank_dir = os.path.join(Configs.outdir, 'blank')
    if not os.path.isdir(blank_dir):
        os.makedirs(blank_dir)

    dirs_to_remove = ['queries', 'sub-alignments']
    if Configs.keep_decomposition:
        dirs_to_remove = []
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
    hmmbuild_paths = []; hmmsearch_paths = []
    if not Configs.hmmdir:
        # default to <outdir>/tree_decomp/root
        Configs.hmmdir = Configs.outdir + '/tree_decomp/root'
    else:
        assert os.path.isdir(Configs.hmmdir), \
                'Provided HMM directory does not exist'

    if not os.path.isdir(Configs.hmmdir):
        # if both backbone alignment/tree are provided by the user
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
            bb_job = BackboneJob(Configs.backbone_path, Configs.query_path,
                    Configs.backbone_tree_path)
            bb_job.setup()

            Configs.backbone_path, Configs.query_path = \
                    bb_job.run_alignment()
            Configs.backbone_tree_path = bb_job.run_tree()

        # after obtaining backbone alignment/tree, perform decomposition
        # and HMMSearches
        print('\nDecomposing the backbone tree...')
        decomp = DecompositionAlgorithm(
                Configs.backbone_path,
                Configs.backbone_tree_path,
                Configs.lower, Configs.alignment_size)
        # legacy version for the first published paper (no differentiation
        # between alignment and assignment levels)
        # NOTE: only run hmmbuilds on sub-alignments in ranges
        if Configs.legacy:
            hmmbuild_paths = decomp.decomposition_legacy(
                    Configs.lower, Configs.upper, lock, pool)
        else:
            hmmbuild_paths = decomp.decomposition(
                    Configs.lower, Configs.upper, lock, pool)
        print('\nPerforming all-against-all HMMSearches ' \
                'between the backbone and queries...')
        search = SearchAlgorithm(hmmbuild_paths)
        hmmsearch_paths = search.search(lock, pool)
        del decomp; del search
    else:
        # go over the given hmm directory and obtain all subset alignment
        # get their retained columns with respect to the backbone alignment
        print('\nFound existing HMM directory: {}'.format(Configs.hmmdir))
        _dummy_search = SearchAlgorithm(None)
        backbone_path = _dummy_search.readHMMDirectory(lock, pool)
        if not Configs.backbone_path:
            Configs.backbone_path = backbone_path

    # create a temp backbone alignment if not availble at <outdir>/tree_decomp/backbone
    tmp_backbone_path, backbone_length = writeTempBackbone(
            Configs.outdir + '/tree_decomp/backbone', Configs.backbone_path)

    ############## codes from WITCH end ###############################

    
    # obtain HMMSearch results
    # hmm_indexes are sorted by their nums of sequences in desending order
    # index to hmms keys=(subset dirname, num seq of the sub-alignment)
    # Optionally using the built paths, which are the ones used
    hmm_indexes, index_to_hmms = obtainHMMs(Configs.hmmdir,
            Configs.lower, Configs.upper, hmmbuild_paths=hmmbuild_paths)

    # obtain weights
    weight_path = Configs.outdir + '/weights.txt'
    if os.path.exists(weight_path):
        print('\nFound existing weights: {}'.format(weight_path))
        scores = readWeightsFromLocal(weight_path)
    else:
        scores = getHMMSearchResults(index_to_hmms)
        # use adjusted bitscore for assignment if specified
        if Configs.use_weight:
            print('\nCalculating weights (adjusted bit-scores)...')
            scores = writeWeights(index_to_hmms, scores, pool) 
        else:
            print('\nLoaded bit-scores...')

        if Configs.save_weight:
            print('\t(user option) Writing weights to {}'.format(weight_path))
            writeWeightsToLocal(scores, weight_path)

    # assign queries to sub-alignments based on ranked bit-scores
    query_assignment = assignQueryToSubset(scores, hmm_indexes, index_to_hmms)

    # write queries to local, split to equal-size subsets if necessary
    query_paths = writeSubQueries(Configs.outdir, hmm_indexes,
            index_to_hmms, Configs.query_path, query_assignment,
            Configs.subproblem_size)

    # run MAFFT-linsi-add for each subproblem
    print('\nSolving each sub-problem with MAFFT-linsi --add...'.format(
        len(query_paths)))
    subalignment_paths = alignSubQueries(Configs.outdir, query_paths)

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
