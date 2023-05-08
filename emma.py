'''
Created on 6.27.2022 by Chengze Shen

Front of eMAFFTadd.
'''

import os, sys, time
from argparse import ArgumentParser, Namespace
#import logging

from src.main import mainAlignmentProcess
from src import * 
from configs import *

version = "0.1.0"

def main():
    parser = _init_parser()
    cmdline_args = sys.argv[1:]

    global main_config_path    
    opts = Namespace()
    main_cmd_defaults = []
    
    # TODO: setup main configuration file
    #if os.path.exists(main_config_path):
    #    with open(main_config_path, 'r') as cfile:
    #        main_cmd_defaults = _read_config_file(cfile, opts)
    input_args = main_cmd_defaults + cmdline_args
    args = parser.parse_args(input_args, namespace=opts)

    buildConfigs(args)
    getConfigs()

    Configs.log('eMAFFTadd is running with: {}'.format(' '.join(cmdline_args)))
    if os.path.exists(main_config_path):
        Configs.log('Main configuration loaded from {}'.format(
            main_config_path))
    
    mainAlignmentProcess(args)


def _init_parser():
    parser = ArgumentParser(description=(
        "This program runs eMAFFTadd, an alignment method "
        "that aligns a portion of sequences first, and "
        "adding in the remaining sequences to the backbone alignment. "
        "The goal is scaling MAFFT-linsi-add, an accurate "
        "sequence adding technique that cannot run on large datasets."),
        conflict_handler='resolve')
    parser.add_argument('-v', '--version', action='version',
            version="%(prog)s " + version)
    
    # basic settings
    basic_group = parser.add_argument_group(
            "Basic parameters".upper(),
            ' '.join(["These are basic fields.",
                "Users can choose to provide the backbone alignment/tree,",
                "or even the path to the eHMM from a previous UPP run.",
                "Otherwise, eMAFFTadd will generate them."]))
    parser.groups = dict()
    parser.groups['basic_group'] = basic_group
    basic_group.add_argument('-i', '--input-path', type=str,
            help='Path to the input unaligned file (all sequences).', required=False)
    basic_group.add_argument('-p', '--hmmdir', type=str,
            help='Path to the HMMs directory generated from UPP', required=False)
    basic_group.add_argument('-b', '--backbone-path', type=str,
            help='Path to the backbone alignment', required=False)
    basic_group.add_argument('-e', '--backbone-tree-path', type=str,
            help='Path to the backbone tree', required=False)
    basic_group.add_argument('-q', '--query-path', type=str,
            help='Path to the queries file that we want to align', required=False)
    basic_group.add_argument('-d', '--outdir', type=str,
            help='Output directory, default: emafftadd_output/', required=False,
            default='emafftadd_output')
    basic_group.add_argument('-o', '--output-path', type=str,
            help='Output file name, default: est.aln.fasta', required=False,
            default='est.aln.fasta')
    basic_group.add_argument('-t', '--num-cpus', type=int,
            help='Number of cpus for multi-processing, default: -1 (all)',
            required=False, default=-1)
    basic_group.add_argument('--continue-run', action='store_const',
            help=' '.join(['Whether to continue from a previous stopped run.',
                    'The signal is the existence of the following folders:',
                    'queries/, sub-backbones/, sub-alignments/, default off.']),
            default=False, const=True, required=False)
    #basic_group.add_argument('--chunksize', type=int,
    #        help='Chunksize for multiprocessing', required=False,
    #        default=1)

    emafftadd_group = parser.add_argument_group(
            "eMAFFTadd parameters".upper(),
            ' '.join(["These are settings for eMAFFTadd parameters.",
                    "Users can decide the number of sequences in sub-alignments",
                    "that will be used by MAFFT-linsi-add to add the query",
                    "sequences (upper and lower bounds).",
                    "Also users can decide the smallest sequence adding problem."]))
    parser.groups['emafftadd_group'] = emafftadd_group
    emafftadd_group.add_argument('--molecule', type=str,
            help='Whether input is amino/dna/rna, default: dna',
            required=False, default='dna', choices=['amino', 'dna', 'rna'])
    emafftadd_group.add_argument('--lower', type=int, default=50,
            help=' '.join(['The lower bound of number of sequences in',
                    'a sub-alignment to use, default: 50']),
            required=False)
    emafftadd_group.add_argument('--upper', type=int, default=100,
            help=' '.join(['The upper bound of number of sequences in',
                    'a sub-alignment to use, default: 100']),
            required=False)
    emafftadd_group.add_argument('--subproblem-size', type=int, default=500,
            help=' '.join(['The total number of sequences (sub-alignment + query)',
                    'for each sequence adding sub-problem.',
                    'If too many queries are assigned to the same sub-alignment,',
                    'they will be randomly broken down to subsets so that each',
                    'smaller problem will not have more than X sequences,',
                    'default: X=500']),
            required=False)

    return parser

if __name__ == "__main__":
    main()
