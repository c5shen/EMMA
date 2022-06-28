from alignment_tools import Alignment
import numpy as np
from collections import defaultdict
import math, sys, os

# find the HMM within range (size = [100, 500])
def obtainHMMs(indir, lower, upper):
    print('In obtainHMMs', flush=True)
    hmmdir = os.path.join(indir, 'upp_align_txt')
    hmms = os.popen('find {} -name hmmbuild.model.* -type f'.format(hmmdir)).read().split('\n')[:-1]

    hmms_in_range = []
    hmm_indexes = []
    index_to_hmms = {}
    index_to_num_seqs = {}
    for hmm in hmms:
        #print('Reading', hmm)
        lines = os.popen('head -n 15 {}'.format(hmm)).read().split('\n')[:-1]
        split_lines = [l.split() for l in lines]
        _dict = {x[0]: ','.join(x[1:]) for x in split_lines}
        num_seq = int(_dict['NSEQ'])

        if num_seq <= upper and num_seq >= lower:
            hmms_in_range.append(hmm)
            hmm_ind = int('/'.join(hmm.split('/')[:-1]).split('/')[-1].split('_')[-1])
            hmm_indexes.append((hmm_ind, num_seq))

            index_to_num_seqs[hmm_ind] = num_seq
            index_to_hmms[hmm_ind] = hmm

    hmm_indexes = sorted(hmm_indexes, key=lambda x: x[1], reverse=True)
    return hmms_in_range, hmm_indexes, index_to_hmms, index_to_num_seqs

def getHMMSearchResults(hmms):
    print('In getHMMSearchResults', flush=True)
    scores = defaultdict(list)
    for hmm in hmms:
        hmmdir = '/'.join(hmm.split('/')[:-1])
        hmm_ind = int(hmmdir.split('/')[-1].split('_')[-1])

        hmmsearches = os.popen('find {} -name hmmsearch.results.* -type f'.format(
            hmmdir)).read().split('\n')[:-1]

        for res in hmmsearches:
            with open(res, 'r') as f:
                this_scores = eval(f.read())
                for taxon, score in this_scores.items():
                    scores[taxon].append((hmm_ind, score[1]))
    for taxon in scores.keys():
        scores[taxon] = sorted(scores[taxon], key=lambda x: x[1], reverse=True)
    return scores

def assignQueryToSubset(scores, hmm_indexes, index_to_hmms):
    print('In assignQueryToSubset', flush=True)
    query_assignment = defaultdict(list)
    assigned_hmms = set()

    for taxon, score in scores.items():
        hmm_ind = score[0][0]
        query_assignment[hmm_ind].append(taxon)
        assigned_hmms.add(index_to_hmms[hmm_ind])

    return query_assignment, list(assigned_hmms)

def writeSubBackbone(outdir, hmms):
    print('In writeSubBackbone', flush=True)
    subbackbone_dir = os.path.join(outdir, 'sub-backbones')
    if not os.path.isdir(subbackbone_dir):
        os.makedirs(subbackbone_dir)
    
    print(len(hmms))
    for hmm in hmms:
        hmmdir = '/'.join(hmm.split('/')[:-1])
        hmm_ind = int(hmmdir.split('/')[-1].split('_')[-1])

        subaln_path = os.popen('find {} -maxdepth 1 -name hmmbuild.input.* -type f'.format(
            hmmdir)).read().split('\n')[0]
        os.system('cp {} {}/subset_{}_backbone.fasta'.format(
            subaln_path, subbackbone_dir, hmm_ind))


def writeSubQueries(outdir, hmm_indexes, index_to_num_seqs,
        query_path, query_assignment):
    print('In writeSubQueries', flush=True)
    query_dir = os.path.join(outdir, 'queries')
    if not os.path.isdir(query_dir):
        os.makedirs(query_dir)

    query = Alignment(); query.read_file_object(query_path)
    total_num_queries = len(query.keys())
    ind_to_query = {}
    for subset_ind, taxa in query_assignment.items():
        taxa_that_exist = set(taxa).intersection(set(query.keys()))
        sub_query = query.sub_alignment(taxa_that_exist)
        ind_to_query[subset_ind] = sub_query
        for t in taxa_that_exist:
            query.pop(t)

    # if any query sequence is not aligned with any subset, put it to the largest
    # subset (arbitrarily tie-break)
    if len(query) > 0:
        for t, s in query.items():
            _hmm_ind = hmm_indexes[0][0]
            if _hmm_ind in ind_to_query:
                ind_to_query[_hmm_ind][t] = s
            else:
                ind_to_query[_hmm_ind] = Alignment()
                ind_to_query[_hmm_ind][t] = s
            #query_assignment[hmm_indexes[0][0]][t] = s
    
    counted_queries = 0
    for _hmm_ind, a in ind_to_query.items():
        counted_queries += len(a.keys())
    assert counted_queries == total_num_queries

    query_paths = []
    for subset_ind, sub_aln in ind_to_query.items():
        # make sure each subset+query combo has <= 500 taxa in total
        query_size = len(sub_aln); bb_size = index_to_num_seqs[subset_ind] 
        if query_size == 0:
            continue
        if bb_size >= 500:
            cutoff = bb_size + 200
        else:
            cutoff = 500
        problem_size = query_size + bb_size 
        num_query_sets = math.ceil(query_size / (cutoff - bb_size))
        subproblem_size = math.ceil(query_size / num_query_sets)

        # shuffle the query list and partition it to [num_query_sets] sets
        query_names = list(sub_aln.keys())
        np.random.shuffle(query_names)
        #print(query_names, problem_size, num_query_sets, file=sys.stderr)
        for i in range(0, num_query_sets, 1):
            lower = i * subproblem_size
            upper = min((i+1) * subproblem_size, query_size)
            subproblem_query_names = query_names[lower:upper]
            subproblem_queries = sub_aln.sub_alignment(subproblem_query_names)

            query_path = os.path.join(query_dir, 'subset_{}_query_{}.fasta'.format(
                subset_ind, i))
            subproblem_queries.write(query_path, 'FASTA')
            query_paths.append(query_path)
    print('Total number of sub-problems: {}'.format(len(query_paths)), flush=True)
    return query_paths

def runMafftAdd(outdir, t, query_paths, hmm_indexes):
    print('In runMafftAdd', flush=True)
    binary = '/home/chengze5/tallis/softwares/bin/mafft-linsi'
    add_dir = os.path.join(outdir, 'sub-alignments')
    if not os.path.isdir(add_dir):
        os.makedirs(add_dir)

    out_paths = []
    for query_path in query_paths:
        filename = query_path.split('/')[-1]
        parts = filename.split('.')[0].split('_')
        subset_ind, subproblem_ind = int(parts[1]), int(parts[3])
        aln_path = os.path.join(outdir, 'sub-backbones',
                'subset_{}_backbone.fasta'.format(subset_ind))
        if not (os.path.exists(aln_path) and os.path.exists(query_path)):
            continue
        out_path = os.path.join(add_dir,
                'subset_{}_query_{}.est.aln.fasta'.format(subset_ind,
                    subproblem_ind))
        if not (os.path.exists(out_path) and os.stat(out_path).st_size != 0):
            cmd = '{} --quiet --thread {} --add {} {} > {}'.format(
                    binary, t, query_path, aln_path, out_path)
            os.system(cmd)
        aln = Alignment(); aln.read_file_object(out_path)
        for taxon in aln.keys():
            aln[taxon] = aln[taxon].upper()
        aln.write(out_path, 'FASTA')
        out_paths.append(out_path)
    return out_paths

    #for ind in hmm_indexes:
    #    aln_path = os.path.join(outdir, 'sub-backbones',
    #            'subset_{}_backbone.fasta'.format(ind[0]))
    #    query_path = os.path.join(outdir, 'queries',
    #            'subset_{}_query.fasta'.format(ind[0]))
    #    if not (os.path.exists(aln_path) and os.path.exists(query_path)):
    #        continue
    #    out_path = os.path.join(add_dir,
    #            'subset_{}.est.aln.fasta'.format(ind[0]))
    #    cmd = '{} --quiet --thread {} --add {} {} > {}'.format(
    #            binary, t, query_path, aln_path, out_path)
    #    os.system(cmd)
    #    aln = Alignment(); aln.read_file_object(out_path)
    #    for taxon in aln.keys():
    #        aln[taxon] = aln[taxon].upper()
    #    aln.write(out_path, 'FASTA')

def runMerger(outdir, aln_paths, backbone_path, t):
    print('In runMerger', flush=True)
    # give an order file
    add_dir = os.path.join(outdir, 'sub-alignments')
    #aln_paths = os.popen('ls {}'.format(add_dir)).read().split('\n')[:-1]

    # first copy over the full backbone path
    os.system('cp {} {}/backbone.fasta'.format(backbone_path, add_dir))
    bb_path = 'backbone.fasta'
    order_path = os.path.join(outdir, 'order.txt')
    order_file = open(order_path, 'w')
    order_file.write('{}/{}\n'.format(add_dir, bb_path))

    for p in aln_paths:
        order_file.write('{}\n'.format(p))
    order_file.close()

    # run merger
    binary = '/home/chengze5/tallis/softwares/ATMerger/merger.py'
    cmd = 'python3 {} -d {} -o {}/est.aln.fasta --order {} -t {}'.format(
            binary, add_dir, outdir, order_path, 1)

    os.system(cmd)
    
def main():
    #upp_dir = '/home/chengze5/tallis/aln_with_long_seq/benchmark_results/ehmm/full50/1000M1-indel/R0/full'
    #outdir = 'temp'; query_path = None; backbone_path = None
    upp_dir = sys.argv[1]
    backbone_path = sys.argv[2]
    query_path = sys.argv[3]
    outdir = sys.argv[4]
    lower, upper = int(sys.argv[5]), int(sys.argv[6])
    t = int(sys.argv[7])
    if len(sys.argv) > 8:
        cont = sys.argv[8] == "continue"
    else:
        cont = False

    assert os.path.isdir(upp_dir)
    assert os.path.exists(backbone_path)
    assert os.path.exists(query_path)

    hmms, hmm_indexes, index_to_hmms, index_to_num_seqs \
            = obtainHMMs(upp_dir, lower, upper)
    if not cont:
        scores = getHMMSearchResults(hmms)
        
        query_assignment, assigned_hmms = assignQueryToSubset(scores,
                hmm_indexes, index_to_hmms)

        writeSubBackbone(outdir, assigned_hmms)

        query_paths = writeSubQueries(outdir, hmm_indexes, index_to_num_seqs,
                query_path, query_assignment)
    else:
        query_paths = os.popen('ls {}/queries'.format(outdir)).read().split('\n')[:-1]
        query_paths = [outdir + '/queries/{}'.format(x) for x in query_paths]

    out_paths = runMafftAdd(outdir, t, query_paths, hmm_indexes)

    runMerger(outdir, out_paths, backbone_path, t)


if __name__ == "__main__":
    main()
