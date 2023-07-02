from configs import Configs
from helpers.alignment_tools import Alignment
import time, os

'''
Align each sub-problem by adding assigned query sequences to the target
sub-alignment using MAFFT-linsi-add
'''
def alignSubQueries(outdir, query_paths):
    Configs.log('Aligning each sub-problem of adding query sequences ' + \
            'to the target sub-alignment')
    start = time.time()

    aln_dir = os.path.join(outdir, 'sub-alignments')
    if not os.path.isdir(aln_dir):
        os.makedirs(aln_dir)

    out_paths = []
    for path in query_paths:
        # from filename to deduce the sub-alignment to use
        filename = path.split('/')[-1]
        parts = filename.split('.')[0].split('_')
        alignment_ind, subproblem_ind = int(parts[1]), int(parts[3])
        #subaln_path = os.path.join(outdir, 'sub-backbones',
        #        'subset_{}_backbone.fasta'.format(hmm_ind))
        subaln_path = os.path.join(outdir, 'tree_decomp/root',
                'A_{}'.format(alignment_ind),
                'subset.aln.fasta')
        if not (os.path.exists(subaln_path) and os.path.exists(path)):
            raise ValueError('Either {} or {} not found!'.format(subaln_path,
                path))

        out_path = os.path.join(aln_dir, 'subset_{}_query_{}.est.aln.fasta'.format(
            alignment_ind, subproblem_ind))
        if not (os.path.exists(out_path) and os.stat(out_path).st_size != 0):
            cmd = '{} --localpair --maxiterate 1000 --quiet --thread {} --add {} {} > {}'.format(
                    Configs.mafftpath, Configs.num_cpus, path, subaln_path, out_path)
            Configs.debug('[MAFFT-add] Command used: {}'.format(cmd))
            os.system(cmd)

        # ensure outputs are all upper cases
        aln = Alignment(); aln.read_file_object(out_path)
        for taxon in aln.keys():
            aln[taxon] = aln[taxon].upper()
        aln.write(out_path, 'FASTA')
        del aln
        out_paths.append(out_path)
        msg = 'Finished sub-problem {}-{}, output: {}'.format(
            alignment_ind, subproblem_ind, out_path)
        Configs.log(msg); print(msg)

    time_aln = time.time() - start
    Configs.log('Done aligning all sub-problems.')
    Configs.runtime('Time to align all sub-problems (s): {}'.format(time_aln))
    return out_paths
