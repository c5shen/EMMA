'''
Created on 6.27.2022 by Chengze Shen

Configuration of eMAFFTadd.
'''

import os, time
try:
    import configparser
except ImportError:
    import ConfigParser as configparser
from src import *
from argparse import ArgumentParser, Namespace
from platform import platform
import shutil

_root_dir = os.path.dirname(os.path.realpath(__file__))
main_config_path = os.path.join(_root_dir, 'main.config')

'''
Configurations defined by users
'''
class Configs:
    global _root_dir

    hmmdir = None
    input_path = None
    backbone_path = None
    backbone_tree_path = None
    query_path = None
    outdir = None
    output_path = None

    lower = 50
    upper = 100
    subproblem_size = 500
    num_cpus = -1

    # hmmalign/hmmsearch/magus paths
    magus_path = os.path.join(_root_dir, 'tools/magus/magus.py')
    if 'macOS' in platform():
        hmmer_dir = os.path.join(_root_dir, 'tools/macOS')
        fasttreepath = os.path.join(_root_dir, 'tools/macOS/FastTreeMP')
        mafftpath = os.path.join(_root_dir, 'tools/macOS/mafft')
    else:
        hmmer_dir = os.path.join(_root_dir, 'tools/hmmer')
        fasttreepath = os.path.join(_root_dir, 'tools/fasttree/FastTreeMP')
        mafftpath = os.path.join(_root_dir, 'tools/mafft')
    hmmalignpath = os.path.join(hmmer_dir, 'hmmalign')
    hmmsearchpath = os.path.join(hmmer_dir, 'hmmsearch')
    hmmbuildpath = os.path.join(hmmer_dir, 'hmmbuild')
    
    ## mafft binary - use the one users installed if possible
    #if shutil.which('mafft'):
    #    mafftpath = shutil.which('mafft')


    log_path = None
    error_path = None
    debug_path = None
    runtime_path = None

    @staticmethod
    def warning(msg, path=None):
        path = Configs.log_path if path is None else path
        Configs.write(msg, 'WARNING', path)

    @staticmethod
    def log(msg, path=None):
        path = Configs.log_path if path is None else path
        Configs.write(msg, 'LOG', path)

    @staticmethod
    def debug(msg, path=None):
        path = Configs.debug_path if path is None else path
        Configs.write(msg, 'DEBUG', path)

    @staticmethod
    def error(msg, path=None):
        path = Configs.error_path if path is None else path
        Configs.write(msg, 'ERROR', path)
    
    @staticmethod
    def runtime(msg, path=None):
        path = Configs.runtime_path if path is None else path
        with open(path, 'a') as f:
            f.write('{}\n'.format(msg))

    @staticmethod
    def write(msg, level, path):
        if path is not None:
            with open(path, 'a') as f:
                f.write('{}\t[{}] {}\n'.format(time.strftime('%Y-%m-%d %H:%M:%S'),
                    level, msg))

# print a list of all configurations
def getConfigs():
    print('\n********* Configurations **********')
    for k, v in Configs.__dict__.items():
        if valid_attribute(k, v):
            print('\tConfigs.{}: {}'.format(k, v))

# valid attribute check
def valid_attribute(k, v):
    assert isinstance(k, str)
    if isinstance(v, staticmethod):
        return False
    if not k.startswith('_'):
        return True
    return False

'''
Build configurations
'''
def buildConfigs(args):
    Configs.input_path = os.path.realpath(args.input_path)
    Configs.hmmdir = os.path.realpath(args.hmmdir)
    Configs.backbone_path = os.path.realpath(args.backbone_path)
    Configs.backbone_tree_path = os.path.realpath(args.backbone_tree_path)
    Configs.query_path = os.path.realpath(args.query_path)
    
    Configs.outdir = os.path.realpath(args.outdir)
    if not os.path.exists(Configs.outdir):
        os.makedirs(Configs.outdir)
    Configs.output_path = os.path.join(Configs.outdir, args.output_path)

    Configs.log_path = os.path.join(Configs.outdir, 'log.txt')
    Configs.error_path = os.path.join(Configs.outdir, 'error.txt')
    Configs.debug_path = os.path.join(Configs.outdir, 'debug.txt')
    Configs.runtime_path = os.path.join(Configs.outdir, 'runtime_breakdown.txt')

    #Configs.chunksize = args.chunksize
    if args.num_cpus > 0:
        Configs.num_cpus = args.num_cpus
    else:
        Configs.num_cpus = os.cpu_count()

    # emafftadd settings
    Configs.molecule = args.molecule
    Configs.lower = args.lower
    Configs.upper = args.upper
    Configs.subproblem_size = args.subproblem_size

    # add any additional arguments to Configs
    #for k in args.__dict__.keys():
    #    if k not in Configs.__dict__:
    #        k_attr = getattr(args, k)

    #        # check whether the configuration is valid
    #        set_valid_configuration(k, k_attr)
