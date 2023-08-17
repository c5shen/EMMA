'''
Created on 1.19.2022 by Chengze Shen

Allow backbone alignment (by MAGUS) and backbone tree estimation (by FastTree2)
'''

import os, subprocess, time
import random
from src import * 
from configs import Configs, valid_attribute
from helpers.alignment_tools import MutableAlignment
from argparse import Namespace

'''
Backbone alignment/tree job
'''
class BackboneJob(object):
    # default setting for a backbone alignment job
    def __init__(self, backbone_path=None, query_path=None,
            backbone_tree_path=None):
        self.alignment_method = 'magus'
        self.backbone_size = None 
        self.backbone_threshold = 0.25
        self.selection_strategy = 'median_length' 
        self.alignment_path = Configs.magus_path

        self.tree_method = 'FastTree2'
        self.tree_path = Configs.fasttreepath

        self.unaligned_backbone_path = None
        self.backbone_path = backbone_path
        self.query_path = query_path
        self.backbone_tree_path = backbone_tree_path

        # magus options
        self.magus_options = Namespace()

        self.outdir = Configs.outdir + '/tree_decomp/backbone'

    # set up from Configs.backbone namespace
    def setup(self):
        if getattr(Configs, 'Backbone', None) != None:
            for k, v in Configs.Backbone.__dict__.items():
                if v:
                    setattr(self, k, v)
        if getattr(Configs, 'MAGUS', None) != None:
            for k, v in Configs.MAGUS.__dict__.items():
                if v:
                    setattr(self.magus_options, k, v)

        # default to <outdir>/tree_decomp/backbone
        if not os.path.isdir(self.outdir):
            os.makedirs(self.outdir)

        # set backbone and queries paths to default paths
        # (and if the default path has an alignment, use it)
        if not self.unaligned_backbone_path:
            self.unaligned_backbone_path = self.outdir + '/backbone.unaln.fasta'
        if not self.backbone_path:
            self.backbone_path = self.outdir + '/backbone.aln.fasta'
        if not self.query_path:
            self.query_path = self.outdir + '/queries.fasta'
        if not self.backbone_tree_path:
            self.backbone_tree_path = self.outdir + '/backbone.tre'
        
        print('\nUsing the following settings for the backbone:')
        for k, v in self.__dict__.items():
            if valid_attribute(k, v):
                if k == 'backbone_size' and v is None:
                    print('\tBackboneJob.{}: min(1000, len(taxa))'.format(k))
                else:
                    print('\tBackboneJob.{}: {}'.format(k, v))

    # split sequences to backbone/query based on the selection strategy
    def splitSequences(self, sequences):
        assert isinstance(sequences, MutableAlignment)

        seq_lengths = sorted([len(seq) for seq in sequences.values()])
        lengths = len(seq_lengths)

        if self.backbone_size is None:
            self.backbone_size = min(1000, lengths)
        else:
            self.backbone_size = int(self.backbone_size)
        Configs.log('Backbone size set to: {}'.format(self.backbone_size))
        backbone_sequences, queries = MutableAlignment(), MutableAlignment()

        if self.selection_strategy == 'median_length':
            l2 = int(lengths / 2)
            if lengths % 2 == 1 or l2 == lengths - 1:
                median_full_length = seq_lengths[l2]
            else:
                median_full_length = (seq_lengths[l2] + seq_lengths[l2+1]) / 2.0
            min_length = int(median_full_length * (1 - self.backbone_threshold))
            max_length = int(median_full_length * (1 + self.backbone_threshold))

            Configs.log('Full length sequences set to be from '
                    + '{} to {} character long'.format(min_length, max_length))

            query_names = [name for name in sequences
                    if len(sequences[name]) > max_length or 
                    len(sequences[name]) < min_length]

            if len(query_names) > 0:
                Configs.log(
                    'Detected {} sequences not within median length'.format(
                        len(query_names)))
                queries = sequences.get_hard_sub_alignment(query_names)
                [sequences.pop(i) for i in list(queries.keys())]

            if len(sequences) < self.backbone_size:
                self.backbone_size = len(sequences)
                Configs.log('Backbone resized to: {}'.format(self.backbone_size))

            sample = sorted(random.sample(
                sorted(list(sequences.keys())), self.backbone_size))
            backbone_sequences = sequences.get_hard_sub_alignment(sample)
            [sequences.pop(i) for i in list(backbone_sequences.keys())]
        elif self.selection_strategy == 'random':
            sample = sorted(random.sample(
                sorted(list(sequences.keys())), self.backbone_size))
            backbone_sequences = sequences.get_hard_sub_alignment(sample)
            [sequences.pop(i) for i in list(backbone_sequences.keys())]
        else:
            Configs.error('Unsupported selection strategy: {}'.format(
                self.selection_strategy))
            notifyError('gcmm/backbone.py')

        # at this point, the asssumption is that backbone_sequences are
        # selected, and some sequences are put in queries
        # WRITE - backbone_sequences -> local
        # MERGE - sequences + queries -> WRITE to local
        unaligned_backbone_path = self.unaligned_backbone_path
        backbone_sequences.write(unaligned_backbone_path, 'FASTA')

        query_path = self.query_path
        sequences.set_alignment(queries)
        sequences.write(query_path, 'FASTA')

        return unaligned_backbone_path, query_path, sequences 
        
    # run alignment
    def run_alignment(self):
        start = time.time()

        if self.alignment_method == 'magus':
            self.path = Configs.magus_path
        
        # some scenarios for interrupted runs
        # (1) backbone/queries are split, and backbone is aligned
        # (2) as (1) but backbone is unaligned
        # (3) backbone/queries are not split
        if os.path.exists(self.backbone_path) \
                and os.stat(self.backbone_path).st_size > 0 \
                and os.path.exists(self.query_path):
            Configs.log('Found existing backbone alignment: {}'.format(
                self.backbone_path))
            print('\nFound existing backbone alignment: {}'.format(
                self.backbone_path))
            # (1.1) query_path is empty meaning all sequences are in backbone
            #       we assume that the user did not stop the run when writing
            #       the query sequences
            if os.stat(self.query_path).st_size == 0:
                print('\nNo query sequences to align! Exiting...')
                Configs.warning('No query sequences to align. Final alignment '
                        + 'saved to {}'.format(Configs.output_path))
                os.system('cp {} {}'.format(self.backbone_path, Configs.output_path))
                exit(0)
            else:
                return self.backbone_path, self.query_path
        elif os.path.exists(self.unaligned_backbone_path) \
                and os.path.exists(self.query_path):
            pass
        else:
            # before we split the backbone/query sequences, make sure
            # the unaligned sequences (input) exist
            if not (Configs.input_path and os.path.exists(Configs.input_path)):
                Configs.error('Input sequences file does not exist')
                raise FileNotFoundError('Input sequences file does not exist')
                #notifyError('gcmm/backbone.py - BackboneJob.run_alignment()')

            # select backbone sequences
            input_sequences = MutableAlignment()
            input_sequences.read_file_object(Configs.input_path)

            self.unaligned_backbone_path, self.query_path, queries = \
                    self.splitSequences(input_sequences)

        # scenarios (2) and (3) both need to run the alignment method,
        # default output dir to <outdir>/tree_decomp/backbone/[method]_alignment
        alignment_outdir = self.outdir + '/{}_alignment'.format(
                self.alignment_method)
        logfile_name = self.outdir \
                + '/{}_alignment_log.txt'.format(self.alignment_method)
        logfile = open(logfile_name, 'w')
        stderrdata, stdoutdata = logfile, logfile

        if self.alignment_method == 'magus':
            cmd = ['python3', self.alignment_path, '--recurse', 'false',
                    '-np', str(Configs.num_cpus),
                    '-i', self.unaligned_backbone_path,
                    '-d', alignment_outdir, '-o', self.backbone_path]
            # load in any presets for MAGUS
            for k, v in self.magus_options.__dict__.items():
                if v:
                    cmd.extend(['--{}'.format(k), str(v)])
        elif self.alignment_method == 'mafft':
            cmd = [self.alignment_path, '--quiet',
                    '--thread', str(Configs.num_cpus),
                    self.unaligned_backbone_path]
            stdoutdata = open(self.backbone_path, 'w')

        print('\nRunning {} backbone alignment...'.format(self.alignment_method))
        Configs.log('Running {} backbone alignment...'.format(
            self.alignment_method))
        Configs.debug('[{}] Command used: {}'.format(
            self.alignment_method.upper(), ' '.join(cmd)))
        p = subprocess.Popen(cmd, stdout=stdoutdata, stderr=stderrdata)
        p.wait()
        if not stdoutdata.closed:
            stdoutdata.close()
        if not stderrdata.closed:
            stderrdata.close()
        if not logfile.closed:
            logfile.close()

        # check if backbone alignment is successfully generated
        if not os.path.exists(self.backbone_path):
            Configs.error('Failed to generate {} backbone alignment, '.format(
                self.alignment_method) + 'please check log at {}'.format(
                    logfile_name))
            notifyError('gcmm/backbone.py - BackboneJob.run_alignment()')

        # need to cast all character to upper for MAFFT
        if self.alignment_method == 'mafft':
            a = MutableAlignment(); a.read_file_object(self.backbone_path)
            for taxon in a.keys():
                a[taxon] = a[taxon].upper()
            a.write(self.backbone_path, 'FASTA')

        Configs.log('Finished {} backbone alignment, backbone file: {}'.format(
            self.alignment_method, self.backbone_path) +
            ', query file: {}'.format(self.query_path))
        print('\tBackbone alignment generated at {}'.format(self.backbone_path))
        print('\tQuery sequences saved at {}'.format(self.query_path))

        # PRE-MATURE END - no queries to align. All sequences aligned in
        # the backbone
        if os.stat(self.query_path).st_size == 0:
            print('\nNo query sequences to align! Exiting...')
            Configs.warning('No query sequences to align. Final alignment '
                    + 'saved to {}'.format(Configs.output_path))
            os.system('cp {} {}'.format(self.backbone_path, Configs.output_path))
            exit(0)

        Configs.runtime('Time to align the backbone (s): {}'.format(
            time.time() - start))
        return self.backbone_path, self.query_path

    # run tree estimation
    def run_tree(self):
        start = time.time()

        # some scenarios for backbone tree estimation
        # (1) backbone_tree_path exists (finished)
        # (2) backbone_tree_path does not exist
        if os.path.exists(self.backbone_tree_path) \
                and os.stat(self.backbone_tree_path).st_size > 0:
            Configs.log('Found existing backbone tree: {}'.format(
                self.backbone_tree_path))
            print('\nFound existing backbone tree: {}'.format(
                self.backbone_tree_path))
        else:
            # before estimating the tree, making sure the backbone alignment
            # exists
            if not os.path.exists(self.backbone_path):
                Configs.error('Did not find a backbone alignment when '
                        + 'estimating the backbone tree.')
                notifyError('gcmm/backbone.py - BackboneJob.run_tree()')

            # MP-version of FastTree2
            logfile_name = self.outdir \
                    + '/{}_tree_log.txt'.format(self.tree_method)
            stderrdata = open(logfile_name, 'w')
            stdoutdata = open(self.backbone_tree_path, 'w')

            if Configs.molecule is None:
                Configs.molecule = Configs.inferDataType(self.backbone_path)

            cmd = [self.tree_path, '-gtr']
            if Configs.molecule == 'dna':
                cmd.extend(['-nt'])
            cmd.extend([self.backbone_path])
            print('\nRunning {} backbone tree estimation...'.format(
                self.tree_method))
            Configs.log('Running {} backbone tree estimation...'.format(
                self.tree_method))
            Configs.debug('[{}] Command used: {}'.format(
                self.tree_method.upper(), ' '.join(cmd)))

            os.system('export OMP_NUM_THREADS={}'.format(Configs.num_cpus))
            p = subprocess.Popen(cmd, stdout=stdoutdata, stderr=stderrdata)
            p.wait()
            if not stderrdata.closed:
                stderrdata.close()
            if not stdoutdata.closed:
                stdoutdata.close()
            
            # check if backbone tree is successfully generated
            if not os.path.exists(self.backbone_tree_path):
                Configs.error('Failed to generate {} backbone tree, '.format(
                    self.tree_method) + 'please check log at {}'.format(
                        logfile_name))
                notifyError('gcmm/backbone.py - BackboneJob.run_tree()')

            Configs.log('Finished {} backbone tree, tree file: {}'.format(
                self.tree_method, self.backbone_tree_path))
            print('\tBackbone tree estimation generated at {}'.format(
                self.backbone_tree_path))

            Configs.runtime('Time to estimate the backbone tree (s): {}'.format(
                time.time() - start))

        return self.backbone_tree_path

