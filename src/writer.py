from configs import Configs
import time, os, math, random
from helpers.alignment_tools import Alignment

'''
Write assigned sub-alignments to local (copying from their decomp)
'''
def writeSubAlignment(outdir, assigned_hmms):
    Configs.log('Writing {} assigned sub-alignments to local'.format(
        len(assigned_hmms)))
    start = time.time()

    subaln_dir = os.path.join(outdir, 'sub-backbones')
    if not os.path.isdir(subaln_dir):
        os.makedirs(subaln_dir)

    for item in assigned_hmms:
        hmm_ind, hmmdir = item
        subaln_path = os.popen('find {} -maxdepth 1 -name hmmbuild.input.* -type f'.format(
            hmmdir)).read().split('\n')[0]
        os.system('cp {} {}/subset_{}_backbone.fasta'.format(
            subaln_path, subaln_dir, hmm_ind))

    time_write_subaln = time.time() - start
    Configs.log('Done writing {} assigned sub-alignments to local'.format(
        len(assigned_hmms)))
    Configs.runtime('Time to write(copy) assigned sub-alignments (s): {}'.format(
        time_write_subaln))

'''
Write query sequences to local with respect to their assigned sub-alignment,
split to subsets if necessary
'''
def writeSubQueries(outdir, hmm_indexes, index_to_hmms, query_path,
        query_assignment, subproblem_size):
    Configs.log('Splitting and writing query sequences to local')
    start = time.time()

    query_dir = os.path.join(outdir, 'queries')
    if not os.path.isdir(query_dir):
        os.makedirs(query_dir)

    query = Alignment()
    query.read_file_object(query_path)
    total_num_queries = len(query.keys())

    # assign queries to corresponding sub-alignments (alignment_ind)
    ind_to_query = {}
    for alignment_ind, taxa in query_assignment.items():
        taxa_that_exist = set(taxa).intersection(set(query.keys()))
        sub_query = query.sub_alignment(taxa_that_exist)
        ind_to_query[alignment_ind] = sub_query
        for t in taxa_that_exist:
            query.pop(t)

    # if any query sequences are not assigned to any sub-alignment, then
    # then put it to the largest sub-alignment (arbitrarily tie-break)
    if len(query) > 0:
        for taxon, seq in query.items():
            alignment_ind = hmm_indexes[0][0]
            if alignment_ind in ind_to_query:
                ind_to_query[alignment_ind][taxon] = seq
            else:
                ind_to_query[alignment_ind] = Alignment()
                ind_to_query[alignment_ind][taxon] = seq
    del query

    counted_queries = 0
    for alignment_ind, aln in ind_to_query.items():
        counted_queries += len(aln.keys())
    assert counted_queries == total_num_queries, \
            ('Somehow not all queries are assigned correctly\t' + \
            'assigned query: {}, total query: {}'.format(counted_queries,
                total_num_queries))

    query_paths = []
    for alignment_ind, aln in ind_to_query.items():
        # make sure each sub-alignment + query combo has <= subproblem_size 
        # sequences
        query_size = len(aln)
        backbone_path = os.path.join(outdir, 'tree_decomp/root/A_{}'.format(
            alignment_ind), 'subset.aln.fasta')
        backbone_size = int(os.popen('wc -l {}'.format(
            backbone_path)).read().split(' ')[0]) / 2

        #backbone_size = index_to_hmms[hmm_ind][1]
        if query_size == 0:
            continue
        # if somehow the backbone sub-alignment has more sequences than
        # the cutoff, raise the cutoff to allow at most 200 queries to be
        # aligned with this sub-alignment
        # NOTICE: this is not an intended feature so it can only happen
        #         if the parameters are tweaked by users
        if backbone_size >= subproblem_size:
            cutoff = backbone_size + 200
        else:
            cutoff = subproblem_size

        current_size = query_size + backbone_size
        num_query_subsets = math.ceil(query_size / (cutoff - backbone_size))
        actual_subproblem_size = math.ceil(query_size / num_query_subsets)

        # shuffle the keys and partition queries to subsets
        query_names = list(aln.keys())
        random.shuffle(query_names)

        for i in range(0, num_query_subsets, 1):
            _l = i * actual_subproblem_size
            _r = min((i+1) * actual_subproblem_size, query_size)
            subproblem_query_names = query_names[_l:_r]
            subaln = aln.sub_alignment(subproblem_query_names)

            subproblem_query_path = os.path.join(query_dir,
                    'subset_{}_query_{}.fasta'.format(alignment_ind, i))
            subaln.write(subproblem_query_path, 'FASTA')
            query_paths.append(subproblem_query_path)

    time_write_queries = time.time() - start
    Configs.log('Done splitting and writing query sequences to local')
    Configs.log('Total number of sub-problems to solve: {}'.format(
        len(query_paths)))
    Configs.runtime('Time to split and write queries to sub-queries (s): {}'.format(
        time_write_queries))
    return query_paths
