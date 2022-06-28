from configs import Configs
import time, os
from helpers.alignment_tools import ExtendedAlignment

'''
Merge all sub-alignments to the backbone alignment through transitivity.
All singleton columns from queries will be collapsed together and written
as lower cases.
'''
def mergeAlignments(output_path, subalignment_paths, backbone_path):
    Configs.log('Merging all sub-alignments with transitivity with ' \
            'singletons from queries collapsed')
    start = time.time()

    masked_outpath = output_path + '.masked'
    assert len(subalignment_paths) > 0, 'No sub-alignments found!'

    # read in the backbone alignment
    full_aln = ExtendedAlignment([])
    full_aln.read_file_object(backbone_path)
    full_aln.from_string_to_bytearray()
    backbone_names = {k: 1 for k in full_aln.keys()}

    # read in each subproblem alignment file and mark them up with insertions
    for subalignment_path in subalignment_paths:
        subaln = ExtendedAlignment([])
        query_names, insertions = subaln.read_queries_alignment(
                backbone_names, subalignment_path)
        
        full_aln.merge_in(subaln, False)
        del subaln
    full_aln.from_bytearray_to_string()

    full_aln.write(output_path, 'FASTA')
    Configs.log('Finished merging all sub-alignments, output file: {}'.format(
        output_path))

    # write a masked version to local as well
    full_aln.remove_insertion_columns()
    full_aln.write(masked_outpath, 'FASTA')
    Configs.log('Masked final alignment written to {}'.format(masked_outpath))

    time_merge = time.time() - start
    Configs.runtime('Time to merge all sub-alignments (s): {}'.format(time_merge))
