'''
Created on 1.22.2022 by Chengze Shen
Modified on 6.27.2022 by Chengze Shen
    - reason: to accomodate eMAFFTadd so that only certain sub-alignments
              are used for HMMbuild (those with [lower, upper] number
              of sequences).

Algorithms for tree decomposition and hmmsearch.
'''

import os, subprocess, math, time, psutil, shutil 
from configs import Configs
from helpers.alignment_tools import Alignment, MutableAlignment 
from src.tree import PhylogeneticTree
import dendropy
from dendropy.datamodel.treemodel import Tree
import tempfile
from functools import partial
import re

from helpers.math_utils import lcm
#from multiprocessing import Queue, Lock

'''
***** ADOPTED from sepp/exhaustive.py - ExhaustiveAlgorithm class *****
Class to perform backbone tree decomposition as does in UPP
'''
class DecompositionAlgorithm(object):
    def __init__(self, backbone_path, backbone_tree_path, lower, alignment_size=50):
        self.symfrac = 0.0
        self.ere = 0.59
        self.informat = 'afa'
        self.molecule = Configs.molecule
        self.path = Configs.hmmbuildpath

        self.strategy = 'centroid'              # default in SEPP/UPP
        self.decomp_strategy = 'hierarchical'   # ensemble of HMMs
        if lower > 10:
            self.assignment_size = lower
        else:
            self.assignment_size = 10           # default in UPP

        self.alignment_size = alignment_size

        self.minsubsetsize = 2
        self.pdistance = 1                      # default in SEPP/UPP
        self.distances = {}
        self.maxDiam = None                     # default in SEPP/UPP

        self.backbone_path = backbone_path
        self.backbone_tree_path = backbone_tree_path

        self.outdir = '{}/tree_decomp'.format(Configs.outdir)
        #self.timeout = 30

    '''
    Read in the backbone alignment and tree
    '''
    def read_alignment_and_tree(self):
        Configs.log('Reading backbone alignment: {}'.format(
            self.backbone_path))
        alignment = Alignment()
        alignment.read_file_object(self.backbone_path)

        Configs.log('Reading backbone tree: {}'.format(
            self.backbone_tree_path))
        tree = PhylogeneticTree(dendropy.Tree.get_from_stream(
            open(self.backbone_tree_path, 'r'),
            schema='newick',
            preserve_underscores=True))
        
        return alignment, tree

    '''
    Do the decomposition on the backbone alignment/tree
    Take in a ProcessPoolExecutor for parallelism

    modification on 6.27.2022 - only builds HMMs of certain sizes defined
    '''
    def decomposition(self, lower, upper, lock, pool):
        start = time.time()
        Configs.log('Started decomposing the backbone to eHMM')
        alignment, tree = self.read_alignment_and_tree()

        ## do not decompose smaller than lower (since it is unnecessary)
        #if lower > self.assignment_size:
        #    Configs.log('Updated assignment subset size from {} to {}'.format(
        #        self.assignment_size, lower))
        #    self.assignment_size = lower

        # restrict alignment subset size to at most the total number of leaves
        # in the backbone
        total_backbone_taxa = tree.leaf_node_names()
        if self.alignment_size > len(total_backbone_taxa):
            Configs.log('Updated alignment subset size from {} to {}'.format(
                self.alignment_size, len(total_backbone_taxa)))
            self.alignment_size = len(total_backbone_taxa)

        assert isinstance(alignment, Alignment)
        assert isinstance(tree, PhylogeneticTree)

        tree.get_tree().resolve_polytomies()
        
        # label edges with numbers so that we can assemble them back
        # at the end
        tree.label_edges()

        ##### 6.4.2023 - for different levels of assignment and alignment
        # Alignment will happen at a "higher" level in the decomposition
        # e.g., when subtree is of size 200-400, we define the corresponding
        # HMMs for alignments
        # Assignment will happen at a "lower" level, e.g., when we further
        # decompose the **alignment subtrees** to subtrees of size 50-100,
        # we denote these subtrees for assignment of queries.
        #
        # Once queries are assigned using these assignment subsets, they will
        # be aligned to the corresponding alignment subset (each assignment
        # subset is an absolute subset of exactly one alignment subset).
        alignment_tree_map = PhylogeneticTree(
                Tree(tree.den_tree)).decompose_tree(
                        self.alignment_size,
                        strategy=self.strategy,
                        minSize=self.alignment_size,
                        tree_map={},
                        decomp_strategy='centroid',
                        pdistance=self.pdistance,
                        distances=self.distances,
                        maxDiam=self.maxDiam)
        assert len(alignment_tree_map) > 0, (
                'Alignment Tree could not be decomposed '
                'given the following settings: '
                'decomp_strategy: {}\nalignment_size: {}\nalignment_size: {}'.format(
                    'centroid', self.alignment_size, self.alignment_size))
        Configs.log('Broken into {} alignment subsets.'.format(
            len(alignment_tree_map))) 

        subset_args = []
        outdirprefix = os.path.join(self.outdir, 'root')
        subalignment_problems = []
        for (aln_key, aln_tree) in alignment_tree_map.items():
            assert isinstance(aln_tree, PhylogeneticTree)
            # further decompose each alignment subproblem into assignment
            # subproblems
            Configs.log('Alignment subset {} has {} leaves'.format(
                aln_key, len(aln_tree.leaf_node_names())))
            subalignment_problems.append(aln_key)

            assignment_tree_map = PhylogeneticTree(
                    Tree(aln_tree.den_tree)).decompose_tree(
                            maxSize=self.assignment_size,
                            strategy=self.strategy,
                            minSize=self.minsubsetsize,
                            tree_map={},
                            decomp_strategy=self.decomp_strategy,
                            pdistance=self.pdistance,
                            distances=self.distances,
                            maxDiam=self.maxDiam)
            assert len(assignment_tree_map) > 0, (
                    'Assignment Tree could not be decomposed '
                    'given the following settings: '
                    'decomp_strategy: {}\nassignment_size: {}\nminsubsetsize: {}'.format(
                        self.decomp_strategy, self.assignment_size, self.minsubsetsize))

            Configs.log('Alignment subset {} has {} assignment subsets of sizes: {}'.format(
                aln_key, len(assignment_tree_map),
                ','.join([str(len(x.leaf_node_names())) for x in assignment_tree_map.values()])))
            
            # create the subproblems locally with the following structure:
            # A_[aln_key]
            #   |_<sub-alignment>
            #   |_A_[aln_key]_0
            #   |_A_[aln_key]_1
            # ...

            # first, get a sub-alignment for the current alignment subproblem
            alignment_taxa = aln_tree.leaf_node_names()
            subalignment = alignment.sub_alignment(alignment_taxa)
            aln_dir = os.path.join(outdirprefix, 'A_{}'.format(aln_key))
            if not os.path.isdir(aln_dir):
                os.makedirs(aln_dir)
            subalignment.write(os.path.join(aln_dir, 'subset.aln.fasta'), 'FASTA')

            # then, create the remaining assignment subproblems
            # within (lower, upper) bounds of num_taxa
            alignment_subset_args = []
            for (assign_key, assign_tree) in assignment_tree_map.items():
                subset_taxa = assign_tree.leaf_node_names()
                num_taxa = len(subset_taxa)
                if num_taxa <= upper and num_taxa >= lower:
                    label = 'A_{}/A_{}_{}'.format(
                            aln_key, aln_key, assign_key)
                    subaln = subalignment.sub_alignment(subset_taxa)
                    #subaln.delete_all_gaps()
                    alignment_subset_args.append((label, subaln))

            # if no assignment subsets are within range for this alignment
            # subset, then use the entire alignment subset as one assignment
            # subset.
            if len(alignment_subset_args) == 0:
                Configs.warning('Alignment subset {}'.format(aln_key) + \
                        'cannot be decomposed with desired range: ({}, {})'.format(
                            lower, upper) + \
                        ', using the entire alignment subset as one assignment ' + \
                        'subset.')
                label = 'A_{}/A_{}_0'.format(aln_key, aln_key)
                alignment_subset_args.append((label, subalignment))
            else:
                #Configs.log('Alignment subset {}:'.format(aln_key) + \
                #        ' creating {} subsets of sizes between [{}, {}]'.format(
                #            len(alignment_subset_args), lower, upper))
                Configs.log('Alignment subset {}: '.format(aln_key) + \
                        'creating {} assignment subsets: {}.'.format(
                            len(alignment_subset_args),
                            ','.join([str((x[0], len(x[1]))) for x in alignment_subset_args])))
            subset_args.extend(alignment_subset_args)

        # create all subset alignments and HMMBuild them
        func = partial(subset_alignment_and_hmmbuild, lock, 
                self.path, outdirprefix,
                self.molecule, self.ere, self.symfrac,
                self.informat)
        hmmbuild_paths = list(pool.map(func, subset_args))
        assert len(hmmbuild_paths) == len(subset_args), \
                'Number of HMMs created does not match ' \
                'the original number of subsets'
        Configs.log('Finished creating {} HMMs to {}'.format(
            len(hmmbuild_paths), outdirprefix))
        
        dur = time.time() - start
        Configs.runtime('Time to decompose the backbone (s): {}'.format(
            dur))
        return subalignment_problems, hmmbuild_paths

'''
Class to perform HMMSearch on all hmmbuild subsets and fragment sequences
(Also will break fragments to fragment chunks for better parallelism)
'''
class SearchAlgorithm(object):
    def __init__(self, hmmbuild_paths):
        self.hmmbuild_paths = hmmbuild_paths
        self.unaligned = None
        self.path = Configs.hmmsearchpath

        self.max_chunk_size = 20000             # default in SEPP/UPP

        self.filters = False                    # No filters
        self.elim = 99999999                    # elim value for hmmsearch
        self.piped = False                      # default in SEPP/UPP
        
        self.outdir = Configs.outdir + '/tree_decomp'

    def search(self, lock, pool):
        Configs.log('Running all-against-all HMM searches between queries ' \
                'and HMMs...')
        start = time.time()
        # create fragment chunks and run HMMSearch on all subsets
        # **** USING SEPP exhaustive implementation to find the best number
        # **** of chunks (Least Common Multiple of #subsets & #cpus // #subsets)
        num_chunks = lcm(
                len(self.hmmbuild_paths),
                Configs.num_cpus) // len(self.hmmbuild_paths)
        #num_chunks = Configs.num_cpus
        frag_chunk_paths = self.read_and_divide_unaligned(num_chunks)

        # create subproblems of HMMSearch
        subset_args = []
        for hmmbuild_path in self.hmmbuild_paths:
            hmmsearch_outdir = '/'.join(hmmbuild_path.split('/')[:-1])
            for i in range(0, len(frag_chunk_paths)):
                frag_chunk_path = frag_chunk_paths[i]
                subset_args.append((hmmsearch_outdir, hmmbuild_path,
                    frag_chunk_path, i))

        func = partial(subset_frag_chunk_hmmsearch, lock, self.path,
                self.piped, self.elim, self.filters)
        hmmsearch_paths = list(pool.map(func, subset_args))
        assert len(hmmsearch_paths) == len(subset_args), \
                'It seems that some HMMSearch jobs failed'
        Configs.log('Finished {} HMMSearch jobs.'.format(len(hmmsearch_paths)))

        dur = time.time() - start
        Configs.runtime('Time to run all-against-all HMMSearches (s): {}'.format(
            dur))
        return frag_chunk_paths 

    def read_and_divide_unaligned(self, num_chunks, extra_frags={}):
        max_chunk_size = self.max_chunk_size
        Configs.log('Reading in fragment files and breaking ' + 
                'them into at most {} chunks '.format(num_chunks) +
                '(each having at most {} sequences)'.format(max_chunk_size))

        # read in query paths (either provided or produced by separating
        # original input sequences to backbone and query sequences)
        unaligned = MutableAlignment()
        unaligned.read_file_object(Configs.query_path)
        ids_unaligned = unaligned.keys()

        # test if input fragment names contain whitespaces or tabs
        ids_unaligned_violations = [id_ for id_ in ids_unaligned
                if (' ' in id_) or ('\t' in id_)]
        if len(ids_unaligned_violations) > 0:
            raise ValueError(
                "Your input fragment file contains {} sequences, ".format(
                    len(ids_unaligned_violations)) + 
                "which names contain either whitespaces or tabs '\\t'. " +
                "Their names are:\n {}".format(
                    "'\n'  ".join(ids_unaligned_violations)))

        # add in extra fragments if provided
        for k, v in extra_frags.items():
            unaligned[k] = v.replace("-", "")

        alg_chunks = unaligned.divide_to_equal_chunks(num_chunks,
                max_chunk_size)
        
        # write to local (<WITCH_OUTPUT_DIR>/tree_decomp/fragment_chunks)
        # and save all frag chunk paths
        ret = []
        fc_outdir = self.outdir + '/fragment_chunks'
        if not os.path.isdir(fc_outdir):
            os.makedirs(fc_outdir)

        for i in range(0, len(alg_chunks)):
            temp_file = None
            if alg_chunks[i]:
                temp_file = '{}/fragment_chunk_{}.fasta'.format(fc_outdir, i)
                #temp_file = tempfile.mktemp(
                #        prefix='fragment_chunk_{}'.format(i),
                #        suffix='.fasta', dir=fc_outdir)
                alg_chunks[i].write(temp_file, 'FASTA')
                Configs.debug('Writing alignment chunk #{} to {}'.format(
                    i, temp_file))
                ret.append(temp_file)
        Configs.log("Finished breaking fragments into {} chunks to {}".format(
            len(ret), fc_outdir))
        return ret


################ HELPER FUNCTIONS FOR MULTI-PROCESSING ################

'''
Obtain subset alignment given taxa, and run hmmbuild on the subset
alignment.
'''
def subset_alignment_and_hmmbuild(lock, binary, outdirprefix, molecule, 
        ere, symfrac, informat, args):
    label, subalignment = args
    subalignment.delete_all_gaps()

    outdir = os.path.join(outdirprefix, label)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    # break down label for file names
    label = label.split('/')[-1]
    
    # write subalignment to outdir
    subalignment_path = '{}/hmmbuild.input.{}.fasta'.format(outdir, label)
    #subalignment_path = tempfile.mktemp(prefix='hmmbuild.input.',
    #        suffix='.fasta', dir=outdir)
    subalignment.write(subalignment_path, 'FASTA')

    # run HMMBuild with 1 cpu given the subalignment
    hmmbuild_path = '{}/hmmbuild.model.{}'.format(outdir, label)
    #hmmbuild_path = tempfile.mktemp(prefix='hmmbuild.model.',
    #        dir=outdir)
    cmd = [binary, '--cpu', '1',
            '--{}'.format(molecule),
            '--ere', str(ere),
            '--symfrac', str(symfrac),
            '--informat', informat,
            '-o', '/dev/null']
    cmd.extend([hmmbuild_path, subalignment_path])
    os.system(' '.join(cmd))
    lock.acquire()
    try:
        Configs.debug('[HMMBuild] Command used: {}'.format(' '.join(cmd)))
    finally:
        lock.release()
        return hmmbuild_path


'''
a single HMMSearch job between a frag chunk and an hmm
'''
def subset_frag_chunk_hmmsearch(lock, binary, piped, elim, filters, args):
    outdir, hmm, unaligned, frag_index = args
    hmm_label = hmm.split('/')[-1].split('.')[-1]
    
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    hmmsearch_path = '{}/hmmsearch.results.{}.fragment_chunk_{}'.format(
            outdir, hmm_label, frag_index)
    #hmmsearch_path = tempfile.mktemp(
    #        prefix='hmmsearch.results.fragment_chunk_{}.'.format(frag_index),
    #        dir=outdir)
    cmd = [binary, '--cpu', '1', '--noali', '-E', str(elim)]
    if not piped:
        cmd.extend(['-o', hmmsearch_path])
    if not filters:
        cmd.extend(['--max'])
    cmd.extend([hmm, unaligned])
    os.system(' '.join(cmd))

    # modify the output file and only retain taxon name, E-value, bit-score
    res = evalHMMSearchOutput(hmmsearch_path)
    with open(hmmsearch_path, 'w') as f:
        f.write(str(res))

    lock.acquire()
    try:
        Configs.debug('[HMMSearch] Command used: {}'.format(' '.join(cmd)))
    finally:
        lock.release()
        return hmmsearch_path

'''
helper function for modifying HMMSearch output file
'''
def evalHMMSearchOutput(path):
    outfile = open(path, 'r')
    results = {}

    pattern = re.compile(
        r"([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+"
        r"([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)")
    start_reading = False
    for line in outfile:
        line = line.strip()
        if not start_reading and line.startswith("E-value") is True:
            start_reading = True
        elif start_reading and line == "":
            start_reading = False
            break
        elif start_reading:
            matches = pattern.search(line)
            if matches is not None and matches.group(0).find("--") == -1:
                results[matches.group(9).strip()] = (
                    float(matches.group(1).strip()),
                    float(matches.group(2).strip()))
                # _LOG.debug("Fragment scores;"
                #           "fragment:%s E-Value:%s BitScore:%s" %(matches
                # .group(9).strip(),matches.group(1).strip(), matches.
                # group(2).strip()))
    outfile.close()
    return results
