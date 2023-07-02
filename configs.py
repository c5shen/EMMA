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

    legacy = False
    use_weight = False
    lower = 50
    upper = 100
    alignment_size = 50
    subproblem_size = 500
    num_cpus = -1
    continue_run = False

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
    
    # mafft binary - use the one users installed if possible
    if shutil.which('mafft'):
        mafftpath = shutil.which('mafft')


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

# check for valid configurations and set them
def set_valid_configuration(name, conf):
    assert isinstance(conf, Namespace), \
            'Looking for Namespace object but find {}'.format(type(conf))

    # backbone alignment settings
    if name == 'Backbone':
        for k in conf.__dict__.keys():
            attr = getattr(conf, k)
            if not attr:
                continue

            if k == 'alignment_method':
                assert str(attr).lower() in ['magus', 'mafft'], \
                    'Backbone alignment method {} not implemented'.format(attr)
            elif k == 'backbone_size':
                assert int(attr) > 0, 'Backbone size needs to be > 0'
            elif k == 'selection_strategy':
                assert str(attr).lower() in ['median_length', 'random'], \
                    'Selection strategy {} not implemented'.format(attr)
            elif k == 'path':
                assert os.path.exists(os.path.realpath(str(attr))), \
                    '{} does not exist'.format(os.path.realpath(str(attr)))
        setattr(Configs, name, conf)
    # settings that change basic Configs class variables such as:
    # fasttreepath, hmmalignpath, etc.
    elif name == 'Basic':
        for k in conf.__dict__.keys():
            attr = getattr(conf, k)
            if not attr:
                continue
            # set variable [k] to [attr] if provided
            setattr(Configs, k, attr)
    elif name == 'MAGUS':
        setattr(Configs, name, conf)

# valid attribute check
def valid_attribute(k, v):
    assert isinstance(k, str)
    if isinstance(v, staticmethod):
        return False
    if not k.startswith('_'):
        return True
    return False

# print a list of all configurations
def getConfigs():
    print('\n********* Configurations **********')
    for k, v in Configs.__dict__.items():
        if valid_attribute(k, v):
            print('\tConfigs.{}: {}'.format(k, v))

'''
Read in from config file if it exists. Any cmd-line provided configs will
override the config file.

Original functionality comes from SEPP -> sepp/config.py
'''
def _read_config_file(filename, opts, expand=None):
    Configs.debug('Reading config from {}'.format(filename))
    config_defaults = []
    cparser = configparser.ConfigParser()
    cparser.optionxform = str
    cparser.read_file(filename)

    if cparser.has_section('commandline'):
        for k, v in cparser.items('commandline'):
            config_defaults.append('--{}'.format(k))
            config_defaults.append(v)

    for section in cparser.sections():
        if section == 'commandline':
            continue
        if getattr(opts, section, None):
            section_name_space = getattr(opts, section)
        else:
            section_name_space = Namespace()
        for k, v in cparser.items(section):
            if expand and k == 'path':
                v = os.path.join(expand, v)
            section_name_space.__setattr__(k, v)
        opts.__setattr__(section, section_name_space)
    return config_defaults

'''
Build configurations
'''
def buildConfigs(args):
    if args.input_path != None:
        Configs.input_path = os.path.realpath(args.input_path)
    if args.hmmdir != None:
        Configs.hmmdir = os.path.realpath(args.hmmdir)
    if args.backbone_path != None:
        Configs.backbone_path = os.path.realpath(args.backbone_path)
    if args.backbone_tree_path != None:
        Configs.backbone_tree_path = os.path.realpath(args.backbone_tree_path)
    if args.query_path != None:
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
    Configs.legacy = args.legacy
    Configs.use_weight = args.use_weight != 0
    Configs.molecule = args.molecule
    Configs.lower = args.lower
    Configs.upper = args.upper
    Configs.alignment_size = args.alignment_size
    Configs.subproblem_size = args.subproblem_size

    Configs.continue_run = args.continue_run

    # add any additional arguments to Configs
    for k in args.__dict__.keys():
        if k not in Configs.__dict__:
            k_attr = getattr(args, k)

            # check whether the configuration is valid
            set_valid_configuration(k, k_attr)
            #setattr(Configs, k, k_attr)
