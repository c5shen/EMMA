import time, os
from configs import Configs, tqdm_styles
from helpers.alignment_tools import Alignment

from tqdm import tqdm

'''
Helper function to return a list of finished and unfinished paths
'''
def updateJobs(outdir, aln_dir, finished, unfinished):
    out_paths, queued_paths = [], []

    for item in finished:
        parts = item.split('-')
        out_paths.append('{}/subset_{}_query_{}.est.aln.fasta'.format(
            aln_dir, parts[0], parts[1]))
    for item in unfinished:
        parts = item.split('-')
        queued_paths.append('{}/queries/subset_{}_query_{}.fasta'.format(
            outdir, parts[0], parts[1]))
    return out_paths, queued_paths

'''
Helper function to find existing jobs and unfinished jobs from aln_dir
'''
def findExistingJobs(aln_dir):
    files = os.popen('ls {}'.format(aln_dir)).read().split('\n')[:-1]
    if len(files) == 0:
       return [], []
    
    # check for file names (unfinished starts with "temp_")
    finished, unfinished = [], []
    for file in files:
        parts = file.split('.est.aln.fasta')[0].split('_')
        if file.startswith('temp_'):
            unfinished.append('{}-{}'.format(parts[2], parts[4]))
        else:
            finished.append('{}-{}'.format(parts[1], parts[3]))
    return finished, unfinished

'''
Helper function to create args for each sub-problem, from each query path
'''
def getSubproblemArguments(outdir, path):
    aln_dir = os.path.join(outdir, 'sub-alignments')

    # from filename to deduce the sub-alignment to use
    filename = path.split('/')[-1]
    parts = filename.split('.')[0].split('_')
    alignment_ind, subproblem_ind = int(parts[1]), int(parts[3])
    subaln_path = os.path.join(outdir, 'tree_decomp/root',
            'A_{}'.format(alignment_ind),
            'subset.aln.fasta')
    out_path = os.path.join(aln_dir, 'subset_{}_query_{}.est.aln.fasta'.format(
        alignment_ind, subproblem_ind))
    temp_out_path = os.path.join(aln_dir,
            'temp_subset_{}_query_{}.est.aln.fasta'.format(
        alignment_ind, subproblem_ind))

    return alignment_ind, subproblem_ind, subaln_path, out_path, temp_out_path

'''
Align each sub-problem by adding assigned query sequences to the target
sub-alignment using MAFFT-linsi-add
'''
def alignSubQueries(outdir, query_paths):
    Configs.log('Aligning each sub-problem of adding query sequences ' + \
            'to the target sub-alignment')
    start = time.time()

    aln_dir = os.path.join(outdir, 'sub-alignments')
    out_paths, queued_paths = [], query_paths
    if not os.path.isdir(aln_dir):
        os.makedirs(aln_dir)
        # create temp output paths for all sub-problems
        for path in queued_paths:
            _, _, subaln_path, out_path, temp_out_path = getSubproblemArguments(
                    outdir, path)
            os.system('touch {}'.format(temp_out_path)) 
    else:
        # trying to identify whether there are existing output files and
        # unfinished jobs
        print('\nFound existing intermediate solutions to sub-problems: {}'.format(
            aln_dir))
        finished, unfinished = findExistingJobs(aln_dir)
        if len(finished) == 0 and len(unfinished) == 0:
            print('\tEmpty intermediate directory, continue as normal.')
        else:
            print('\tFinished sub-problems ({}): {}'.format(
                len(finished), finished))
            print('\tUnfinished sub-problems ({}): {}'.format(
                len(unfinished), unfinished))
            out_paths, queued_paths = updateJobs(
                    outdir, aln_dir, finished, unfinished)

    assert (len(out_paths) + len(queued_paths)) == len(query_paths), \
            ' '.join(['Number of sub-problems is different between',
                'current run and the jobs found at {}'.format(aln_dir)])

    with tqdm(total=len(queued_paths), **tqdm_styles) as pbar:
        for path in queued_paths:
            alignment_ind, subproblem_ind, subaln_path, out_path, temp_out_path \
                    = getSubproblemArguments(outdir, path)

            if not (os.path.exists(subaln_path) and os.path.exists(path)):
                raise FileNotFoundError('Either {} or {} not found!'.format(subaln_path,
                    path))

            if not (os.path.exists(out_path) and os.stat(out_path).st_size != 0):
                cmd = '{} --localpair --maxiterate 1000 --quiet --thread {} --add {} {} > {}; mv {} {}'.format(
                        Configs.mafftpath, Configs.num_cpus, path,
                        subaln_path, temp_out_path,
                        temp_out_path, out_path)
                Configs.debug('[MAFFT-add] Command used: {}'.format(cmd))
                os.system(cmd)
                #os.system('mv {} {}'.format(temp_out_path, out_path))

            # ensure outputs are all upper cases
            aln = Alignment(); aln.read_file_object(out_path)
            for taxon in aln.keys():
                aln[taxon] = aln[taxon].upper()
            aln.write(out_path, 'FASTA')
            del aln
            out_paths.append(out_path)
            msg = 'Finished sub-problem {}-{}, output: {}'.format(
                alignment_ind, subproblem_ind, out_path)
            Configs.log(msg)#; print(msg)
            pbar.update(1)

    time_aln = time.time() - start
    Configs.log('Done aligning all sub-problems.')
    Configs.runtime('Time to align all sub-problems (s): {}'.format(time_aln))
    return out_paths
