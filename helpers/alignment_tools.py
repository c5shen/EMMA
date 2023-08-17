import time
import numpy as np
import re, os, sys
import itertools
from copy import deepcopy
from random import random
from functools import reduce
from dendropy.datamodel.taxonmodel import Taxon
import time
import io
from collections import defaultdict
try:
    from collections.abc import Mapping
except ImportError:
    from collections import Mapping
from operator import add

try:
    filetypes = (io.IOBase, file)
except NameError:
    filetypes = io.IOBase

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
    from io import BytesIO

gappat = re.compile(r'[-?]+')
nogappat = re.compile(r'[^-?]+')
_INDEL = re.compile(r"[-]")
_COACH_INS = re.compile(r"[.]")
_DANGEROUS_NAME_CHARS = re.compile(r"[^a-zA-Z0-9]")
ILLEGAL_CHARS = re.compile(r"[^a-zA-Z?-]")

# helper function, prints out all functionalities in this python file
def printHelper():
    print("1.\tdef readFastSP(path) --\tRead and parse FastSP output." +
            "\n\treturns a map of SPFN, SPFP, and Expansion.")
    print("2.\tdef readHMMSearch(path) --\tRead and parse HMMSearch output." +
            "\n\treturns a map of taxa to tuples: (e-value, bitscore)")
    print("3.\tdef read_fasta --\tRead in a fasta file, used in Alignment()")
    print("4.\tdef write_fasta --\tWrite an Alignment object to local.")
    print("5.\tclass Alignment() --\tAlignment class, for jobs with alignment")

# Function to read in a log file, which has runtime information attached
# at the bottom (using /usr/bin/time -v to generate)
# return: the runtime in minutes
def readRuntime(file_path):
    if os.path.isfile(file_path):
        runtime = 0
        with open(file_path, 'r') as f:
            lines = f.read().split('\n')
            if len(lines) < 23:
                return np.nan
            raw = lines[-24:-1][4].split('):')[1].strip().split(':')
            if len(raw) == 2:
                runtime = int(raw[0]) + float(raw[1])/60
            elif len(raw) == 3:
                runtime = int(raw[0]) * 60 + int(raw[1]) + float(raw[2])/60
            else:
                runtime = np.nan
        return runtime
    else:
        return np.nan

# Function to read in a fastsp file and parse the results
# return: a map of SPFN, SPFP, Expansion score
def readFastSP(file_path):
    ret = {'SPFN': np.nan, 'SPFP': np.nan, 'error': np.nan,
            'Expansion': np.nan}
    if os.path.isfile(file_path):
        f = open(file_path, 'r')
        lines = f.read().split('\n')[:-1]
    else:
        return ret
    if len(lines) > 0:
        ret['SP-Score'] = float(lines[0].split(' ')[1])
        ret['Modeler-Score'] = float(lines[1].split(' ')[1])
        ret['SPFN'] = float(lines[2].split(' ')[1])
        ret['SPFP'] = float(lines[3].split(' ')[1])
        ret['Expansion'] = float(lines[4].split(' ')[1])
        ret['TC'] = float(lines[5].split(' ')[1])
        ret['error'] = (ret['SPFN'] + ret['SPFP']) / 2
    return ret


################ Copied from PASTA-Alignment        ##################
################ CompactAlignment class definition  ##################
def write_compact_to_fasta(alignment, dest):
    """Writes the `alignment` in FASTA format to either a file object or file"""
    file_obj = None
    if isinstance(dest, str):
        file_obj = open(dest, "w")
    else:
        file_obj = dest
    for name in list(alignment.keys()):
        s = alignment.as_string_sequence(name)
        file_obj.write('>{}\n{}\n'.format(name, s))
    if isinstance(dest, str):
        file_obj.close()

def read_fasta(src, remove_gaps=False):
    """generator that returns (name, sequence) tuples from either a FASTA
    formatted file or file object.
    """
    file_obj = None
    if isinstance(src, str):
        try:
            file_obj = open(src, "r")
        except IOError:
            print(("The file `%s` does not exist, exiting gracefully" % src))
    elif isinstance(src, filetypes):
            file_obj = src
    else:
        raise TypeError('FASTA reader cannot recognize the source of %s, %s, %s' % (src,type(src),isinstance(src, filetypes)))
    name = None
    seq_list = list()
    for line_number, i in enumerate(file_obj):
        if i.startswith('>'):
            if name:
                if remove_gaps:
                    yield name, ''.join(seq_list).replace('-', '')
                else:
                    yield name, ''.join(seq_list)
                seq_list = list()
            name = i[1:].strip()
        else:
            #seq = ''.join(i.strip().upper().split())
            seq = ''.join(i.strip().split())
            #if not is_sequence_legal(seq):
            #    raise Exception("Error: illegal characeters in sequence at line %d" % line_number)
            seq_list.append(seq)
    if name:
        if remove_gaps:
            yield name, ''.join(seq_list).replace('-', '')
        else:
            yield name, ''.join(seq_list)
    if isinstance(src, str):
        file_obj.close()

def write_fasta(alignment, dest):
    """Writes the `alignment` in FASTA format to either a file object or file"""
    file_obj = None
    if isinstance(dest, str):
        file_obj = open(dest, "w")
    else:
        file_obj = dest
    for name, seq in list(alignment.items()):
        file_obj.write('>%s\n%s\n' % (name, seq) )
    if isinstance(dest, str):
        file_obj.close()

'''
Infer data type from file
'''
def inferDataType(path):
    sequences = read_fasta(path, remove_gaps=True)
    acg, t, u, total = 0, 0, 0, 0
    for taxon, seq in sequences:
        letters = seq.upper()
        for letter in letters:
            total = total + 1
            
            if letter in ('A', 'C', 'G', 'N'):
                acg = acg + 1
            elif letter == 'T':
                t = t + 1
            elif letter == 'U':
                u = u + 1
    
    if u == 0 and (acg + t)/total > 0.9:
        #print("Found {}% ACGT-N, assuming DNA..".format(int(100*(acg + t)/total)))
        dataType = "dna"
    elif t == 0 and (acg + u)/total > 0.9:
        #print("Found {}% ACGU-N, assuming RNA..".format(int(100*(acg + u)/total)))
        dataType = "rna"
    else:
        #print("Assuming protein..")
        dataType = "amino"
          
    return dataType


class Alignment(dict, object):
    """A simple class that maps taxa names to sequences.
    TODO: switch to dendropy character_matrix
    """
    def __init__(self):
        "creates an empty matrix"
        dict.__init__(self)
        self.datatype = None

    def get_datatype(self):
        return self._datatype

    def set_datatype(self, d):
        if d is None:
            self._datatype = None
        else:
            self._datatype = d.upper()

    datatype = property(get_datatype, set_datatype)

    def get_sequence_names(self):
        "returns a list of sequence names"
        return list(self.keys())

    def get_num_taxa(self):
        "returns the number sequences"
        return len(self.get_sequence_names())

    def read_file_object(self, file_obj, file_format='FASTA'):
        """Augments the matrix by reading the file object.
        If duplicate sequence names are encountered then the old name will be replaced.
        """
        if ( file_format.upper() == 'FASTA' ):
            read_func = read_fasta
        #elif ( file_format.upper() == 'NEXUS' ):
        #    read_func = read_nexus
        #elif ( file_format.upper() == 'PHYLIP' ):
        #    read_func = read_phylip
        #elif ( file_format.upper() == 'COMPACT3' ):
        #    read_func = read_compact3
        else:
            raise NotImplementedError("Unknown file format (%s) is not supported" % file_format)
        for name, seq in read_func(file_obj):
            self[name] = seq

    def write_filepath(self, filename, file_format='FASTA', zipout=False):
        """Writes the sequence data in the specified `file_format` to `filename`"""
        
        file_obj = open_with_intermediates(filename,'w')
        if zipout:
            file_obj.close()
            file_obj = StringIO()
        self.write(file_obj, file_format=file_format)
        if zipout:
            import gzip
            file_obj_gz = gzip.open(filename, "wb", 6)
            file_obj_gz.write(str.encode(file_obj.getvalue()))
            file_obj_gz.close()
        file_obj.close()

    def write(self, file_obj, file_format):
        """Writes the sequence data in the specified `file_format` to `file_obj`"""
        if ( file_format.upper() == 'FASTA' ):
            write_func = write_fasta
        #elif ( file_format.upper() == 'NEXUS' ):
        #    write_func = write_nexus
        #elif ( file_format.upper() == 'PHYLIP' ):
        #    write_func = write_phylip
        #elif ( file_format.upper() == 'COMPACT' ):
        #    write_func = write_compact            
        #elif ( file_format.upper() == 'COMPACT2' ):
        #    write_func = write_compact2       
        #elif ( file_format.upper() == 'COMPACT3' ):
        #    write_func = write_compact3
        else:
            write_func = write_fasta
        write_func(self, file_obj)

    def unaligned(self):
        """
        Returns a new alignment with all gaps and missing sequences removed.
        """
        new_alignment = Alignment()
        new_alignment.datatype = self.datatype
        for name, seq in self.items():
            new_seq = re.sub(_INDEL, '', str(seq))
            if new_seq != '':
                new_alignment[name] = new_seq
        return new_alignment

    def sub_alignment(self, sub_keys):
        "Creates an new alignment with a subset of the taxa."
        new_alignment = Alignment()
        new_alignment.datatype = self.datatype
        for key in sub_keys:
            if key in self:
                new_alignment[key] = self[key]
        return new_alignment

    def is_empty(self):
        return self.__len__() < 1

    def is_aligned(self):
        if self.is_empty():
            raise ValueError("The alignment is empty.\n")
        else:
            v = list(self.values())
            first_seq_len = len(v[0])
            return all([len(i) == first_seq_len for i in v])

    def partition_info(self, base=0):
        return (self.datatype, 1+base, self.sequence_length() + base)

    def sequence_length(self):
        if self.is_aligned():
            return len(list(self.values())[0])

    def max_sequence_length(self):
        return max(len(re.sub(_INDEL, '', v)) for v in list(self.values()))

    def get_all_gap_cols(self):
        all_gaps = list(range(0, self.sequence_length()))
        for seq in self.values():
            all_gaps[:] = [i for i in all_gaps if seq[i] == '-']
        return all_gaps

    def remove_columns(self, indexes):
        for name, seq in self.items():
            self[name] = ''.join((char for idx, char in enumerate(seq)
                if idx not in indexes))

    def delete_all_gaps(self):
        rem = set(self.get_all_gap_cols())
        subset = [x for x in range(0, self.sequence_length()) if x not in rem]
        self.remove_columns(set(rem))
        return subset
    
    def from_bytearray_to_string(self):
        for k,v in self.items():
            self[k] = str(v)

    def from_string_to_bytearray(self):
        for k,v in self.items():
            self[k] = bytearray(v)   
    
    def divide_to_equal_chunks(self, chunks, max_chunk_size=None):
        names = self.get_sequence_names()
        ret = []
        if max_chunk_size and len(names) / chunks > max_chunk_size:
            chunks = len(names) // max_chunk_size + 1
        for i in range(0, chunks):
            subset = names[i:len(names):chunks]
            if subset:
                subset_alg = self.sub_alignment(subset)
            else:
                subset_alg = None
            ret.append(subset_alg)
        return ret

    def mask_gapy_sites(self,minimum_seq_requirement):        
        n = len(list(self.values())[0])
        # The following implements column-based masking. Seems to be more efficient than row based
        masked = []
        allseqs = list(self.values())
        allseqs.sort(key=lambda x: x.count("-"))
        for c in range(0,n):
            r = minimum_seq_requirement
            for seq in allseqs:
                if seq[c] != "-":
                    r -= 1
                if r <= 0:
                    break
            if r > 0:
                masked.append(c)
                 
        #_LOG.debug("%d Columns identified for masking" %len(masked))
        if not masked:
            return
        included = [z for z in reduce(lambda x,y: x+[(x[-1][1]+1,y)],masked,[(-1,-1)]) if z[0]!=z[1]]
        if not included:
            included.append((masked[-1]+1,n))
        if included[-1][1] < n and masked[-1]+1 != n:
            included.append((masked[-1]+1,n))
        for k,seq in self.items():
            tmp = []
            for (i,j) in included:
                tmp.append(seq[i:j])
            self[k] = "".join(tmp)
        nn = len(list(self.values())[0])
        assert (len(masked) == n - nn), "Masking results is not making sense: %d %d %d" %(len(masked), n , nn)
        #_LOG.debug("Masking done. Before masking: %d; After masking: %d; minimum requirement: %d;" %(n,nn,minimum_seq_requirement))

    def merge_in(self, she):
        merge_in(self,she)


class AlignmentSequence:
    def __init__(self):
        self.seq = None
        self.pos = []
    
    def replace(self, match_char, replace_char):
        s = AlignmentSequence()
        s.pos.extend(self.pos)
        s.seq = self.seq.replace(match_char, replace_char)
        return s
    
    def as_string(self, pad_to):
        s=[]
        nxt = 0
        for i,p in enumerate(self.pos):
            if nxt != p:
                s.append("-" * (p-nxt))
            s.append(self.seq[i]) 
            nxt = p+1
        if (pad_to > nxt):
            s.append("-" * (pad_to-nxt))
        return "".join(s)
    
    def __str__(self):
        return self.as_string(0)
    
    def __repr__(self):
        return self.__str__()

class CompactAlignment(dict,object):
    def __init__(self):
        self.colcount = 0
        self.datatype = None
    
    def sub_alignment(self, sub_keys):
        "Creates an new alignment with a subset of the taxa."
        new_alignment = CompactAlignment()
        new_alignment.datatype = self.datatype
        for key in sub_keys:
            if key in self:
                new_alignment[key] = self[key]
        return new_alignment
   
    def is_aligned(self):
        return True
    def sequence_length(self):
        return self.colcount
    def get_num_taxa(self):
        return len(self)
    
    #def unaligned(self):
    #    n = Alignment()
    #    n.datatype = self.datatype
    #    for k,seq in self.items():
    #        n[k] = seq.seq
    #    return n
    
    def iter_column_character_count(self, seqsubset = None):
        if seqsubset is None:
            seqsubset = list(self.keys())
        
        counts = [0] * self.colcount
        for k in seqsubset:
            for pos in self[k].pos:
                counts[pos] += 1
        
        for c in counts:
            yield c
            
    def iter_columns_with_minimum_char_count(self, x, seqsubset = None):
        for i,c in enumerate(self.iter_column_character_count(seqsubset)):
            if c >= x:
                yield i

    def iter_columns_with_maximum_char_count(self, x, seqsubset = None):
        for i,c in enumerate(self.iter_column_character_count(seqsubset)):
            if c <= x:
                yield i
            
                
    def get_insertion_columns(self,shared):
        return set(i for i in self.iter_columns_with_maximum_char_count(0,shared)) 

    def merge_in(self, she):
        '''
        Merges she inside self, assuming we share some common taxa, and the 
        alignment of common taxa is identical across both alignments.
        
        When assumptions are not met, behavior is largely undefined. 
        '''      
        #global _T_ID
        #_T_ID += 1
        #ID = _T_ID    
        #TIMING_LOG.info("transitivitymerge (%d) started" %ID )    
        mykeys = set(self.keys())
        herkeys = set(she.keys())
        #_LOG.debug("Transitive Merge Started. ID:%d - Rows: %d,%d" %(ID,len(mykeys),len(herkeys)))    
        shared = mykeys.intersection(herkeys)
        #_LOG.debug("Shared seq: %d" %(len(shared)))        
        onlyhers = herkeys - shared
        me_ins = self.get_insertion_columns(shared)
        she_ins = she.get_insertion_columns(shared)
        #_LOG.debug("Insertion Columns: %d,%d" %(len(me_ins),len(she_ins)))
        
        memap=[]
        shemap=[]
        
        melen = self.colcount
        shelen =  she.colcount
    
        ime = 0
        ishe = 0
        inew = 0
        while ime < melen or ishe < shelen:
            if ime in me_ins:              
                memap.append(inew)
                ime += 1
            elif ishe in she_ins:
                shemap.append(inew)
                ishe += 1
            else:       
                memap.append(inew)
                shemap.append(inew)
                ishe += 1
                ime += 1
            inew += 1
            
        self.colcount = inew
        
        for seq in self.values():
            seq.pos = [memap[p] for p in seq.pos]
            
        for k in onlyhers:
            she[k].pos = [shemap[p] for p in she[k].pos]
            self[k] = she[k]        
            
        #TIMING_LOG.info("transitivitymerge (%d) finished" %ID )
        #_LOG.debug("Transitive Merge Finished. ID:%d; cols after: %d" %(ID,self.colcount))

    def mask_gapy_sites(self,minimum_seq_requirement):                
        #_LOG.debug("Masking alignment sites with fewer than %d characters from alignment with %d columns" 
        #           %(minimum_seq_requirement,self.colcount))

        masked = set(self.iter_columns_with_maximum_char_count(minimum_seq_requirement-1))
        
        #_LOG.debug("%d Columns identified for masking" %len(masked))
        
        self.mask_sites(masked)
        
    def mask_unaligned_sites(self):
        #_LOG.debug("Masking alignment sites with lower case letters from an alignment with %d sites" 
        #           %(self.colcount))

        masked = set()
        for seq in self.values():
            for c,i in zip(seq.seq,seq.pos):
                if c > 'a' and c < 'z':
                    masked.add(i)
        
        #_LOG.debug("%d Columns identified for masking" %len(masked))
        
        self.mask_sites(masked)     
        
        return self   
        
    def mask_sites(self, masked):
        if not masked:
            return
        
        off = 0
        colmap = []
        for c in range(0,self.colcount):
            if c in masked:
                off += 1
                colmap.append(-1)
            else:
                colmap.append(c - off)
                
        #_LOG.debug("Column index mapping calculated.")
        for seq in self.values():
            # Find ID of char positions that would be masked
            maskind = [i for i,c in enumerate(seq.pos) if c in masked]
            n = len(seq.pos)
            included = [z for z in reduce(lambda x,y: x+[(x[-1][1]+1,y)],maskind,[(-1,-1)]) if z[0]!=z[1]]
            if not maskind:
                included.append((0,n))
            elif (not included or included[-1][1] < n) and maskind[-1]+1 != n:
                included.append((maskind[-1]+1,n))
            tmp = []
            for (i,j) in included:
                tmp.extend(seq.seq[i:j])
            seq.seq = "".join(tmp)
            seq.pos = [colmap[x] for x in (p for p in seq.pos if p not in masked)]
                
        self.colcount -= off
        #_LOG.debug("Masking done. Before masking: %d; After masking: %d;" 
        #           %(self.colcount+off,self.colcount))
        
    def read_filepath(self, filename, file_format='FASTA'):
        """Augments the matrix by reading the filepath.
        If duplicate sequence names are encountered then the old name will be replaced.
        """
        file_obj = open(filename, 'r')
        ret = self.read_file_object(file_obj, file_format=file_format)
        file_obj.close()
        return ret

    def read_file_object(self, file_obj, file_format='FASTA'):
        """Augments the matrix by reading the file object.
        If duplicate sequence names are encountered then the old name will be replaced.
        """
        if ( file_format.upper() == 'FASTA' ):
            read_func = read_fasta        
        #elif ( file_format.upper() == 'COMPACT' ):
        #    read_func = read_compact
        #elif ( file_format.upper() == 'COMPACT3' ):
        #    read_func = read_compact3
        else:
            raise NotImplementedError("Unknown file format (%s) is not supported" % file_format)
        self.colcount = 0
        for name, seq in read_func(file_obj):
            cseq, l = self.get_alignment_seq_object(seq)
            self[name] = cseq
            self.colcount = max(l, self.colcount)
        
    def as_string_sequence(self,name):
        seq = self[name]
        return seq.as_string(self.colcount)
            
    def get_alignment_seq_object(self, seq):
        cseq = AlignmentSequence()
        l = 0
        if isinstance(seq, str):            
            for m in re.finditer(nogappat, seq):
                cseq.pos.extend(range(m.start(),m.end()))
            cseq.seq = re.sub(gappat,"",seq)
            l = len(seq)
        else:
            cseq.seq = seq[0]
            cseq.pos = seq[1]
            l = seq[1][-1] + 1
        return (cseq,l)

    def update_dict_from(self, alignment):
        for k in self.keys():
            alignment[k] = self.as_string_sequence(k)
        alignment.datatype = self.datatype
            
    def update_from_alignment(self, alignment):
        for k,v in alignment.items():
            self[k],l = self.get_alignment_seq_object(v)
        self.colcount = l
        self.datatype = alignment.datatype
            
    def write_filepath(self, filename, file_format='FASTA', zipout=False):
        """Writes the sequence data in the specified `file_format` to `filename`"""
        
        file_obj = open_with_intermediates(filename,'w')
        if zipout:
            file_obj.close()
            file_obj = StringIO()
        self.write(file_obj, file_format=file_format)
        if zipout:
            import gzip
            file_obj_gz = gzip.open(filename, "wb", 6)
            file_obj_gz.write(str.encode(file_obj.getvalue()))
            file_obj_gz.close()
        file_obj.close()

    def write(self, file_obj, file_format):
        """Writes the sequence data in the specified `file_format` to `file_obj`"""
        if ( file_format.upper() == 'FASTA' ):
            write_func = write_compact_to_fasta        
        #elif ( file_format.upper() == 'COMPACT' ):
        #    write_func = write_compact_to_compact
        #elif ( file_format.upper() == 'COMPACT3' ):
        #    write_func = write_compact_to_compact3 
        #elif ( file_format.upper() == 'PHYLIP' ):
        #    write_func = write_compact_to_phylip 
        else:
            write_func = write_compact_to_fasta
        write_func(self, file_obj)

def compact(alg):
    comp = CompactAlignment()
    comp.update_from_alignment(alg)
    return comp

################# Adapted from SEPP/alignment.py ################
class ReadOnlyAlignment(Mapping, object):
    def get_num_taxa(self):
        return len(self.get_sequence_names())
    
    def get_length(self):
        if self.is_empty():
            raise ValueError("The alignment is empty.\n")
        else:
            return len(next(iter(self.values())))

    def get_sequence_names(self):
        return list(self.keys())

    def write(self, file_obj, file_format):
        if file_format.upper() == 'FASTA':
            write_func = write_fasta
        else:
            write_func = write_fasta
        write_func(self, file_obj)

    def is_empty(self):
        return self.__len__() < 1

    def is_aligned(self):
        if self.is_empty():
            raise ValueError("The alignment is empty.\n")
        else:
            l0 = self.get_length()
            return all(len(i) == l0 for i in self.values())

    def is_all_gap(self, pos):
        """Checks to see if a column is all gap column at position x"""
        for seq in self.values():
            if seq[pos] != '-':
                return False
        return True

    def __str__(self):
        return '\n'.join([">%s\n%s" % (k, self[k])
                          for k in sorted(self.keys())])

    def divide_to_equal_chunks(self, chunks, max_chunk_size=None):
        names = self.get_sequence_names()
        ret = []
        if max_chunk_size and len(names) / chunks > max_chunk_size:
            chunks = len(names) // max_chunk_size + 1
        for i in range(0, chunks):
            subset = names[i:len(names):chunks]
            if subset:
                subset_alg = self.get_soft_sub_alignment(subset)
            else:
                subset_alg = None
            ret.append(subset_alg)
        return ret

    def get_soft_sub_alignment(self, sub_key):
        """
        Returns a read-only sub-alignment, which won't consume extra memory,
        since it will not hold a separate copy of the alignment.
        """
        return ReadonlySubalignment(sub_key, self)

class MutableAlignment(dict, ReadOnlyAlignment, object):
    """ An alignment object, that can be modified. This is the class that
    should be used mainly for holding alignments.
    """
    def __init__(self):
        """creates an empty matrix"""
        dict.__init__(self)
        self.datatype = None

    def set_alignment(self, alignment):
        for name, seq in alignment.items():
            self[name] = seq.upper()

    def read_filepath(self, filename, file_format='FASTA'):
        """Augments the matrix by reading the filepath.
        If duplicate sequence names are encountered then the old name will
        be replaced.
        """
        file_obj = open(filename, 'r')
        return self.read_file_object(file_obj, file_format=file_format)

    def read_file_object(self, file_obj, file_format='FASTA'):
        """Augments the matrix by reading the file object.
        If duplicate sequence names are encountered then the old name will
        be replaced.
        """
        if file_format.upper() == 'FASTA':
            read_func = read_fasta
#        elif (file_format.upper() == 'NEXUS'):
#            read_func = read_nexus
#        elif (file_format.upper() == 'PHYLIP'):
#            read_func = read_phylip
        else:
            raise NotImplementedError(
                "Unknown file format (%s) is not supported" % file_format)
        for name, seq in read_func(file_obj):
            self[name] = seq.upper()
        return self

    def add_column(self, pos, char='-'):
        # _LOG.debug("Added a column to reference alignment at position
        # %d" %pos)
        for name, seq in self.items():
            if hasattr(char, "get"):
                c = char.get(name, '-')
            else:
                c = char
            seq = seq[:pos] + c + seq[pos:]
            self[name] = seq

    def remove_column(self, pos):
        for name, seq in self.items():
            seq = seq[:pos] + seq[pos + 1:]
            self[name] = seq

    def remove_columns(self, indexes):
        for name, seq in self.items():
            self[name] = ''.join((
                char for idx, char in enumerate(seq) if idx not in indexes))

    def keep_columns(self, indexes):
        for name, seq in self.items():
            self[name] = ''.join((
                char for idx, char in enumerate(seq) if idx in indexes))

    def get_all_gap_cols(self):
        all_gaps = list(range(0, self.get_length()))
        for seq in self.values():
            all_gaps[:] = [i for i in all_gaps if seq[i] == '-']
        return all_gaps

    def get_all_nongap_cols(self):
        all_gaps = list(range(0, self.get_length()))
        for seq in self.values():
            all_gaps[:] = [i for i in all_gaps if seq[i] == '-']
        return [x for x in range(0, self.get_length()) if x not in all_gaps]

    def delete_all_gap(self):
        """
        Delete all sites that consists of nothing but gaps
        """
        # pdb.set_trace()

        rem = set(self.get_all_gap_cols())
        subset = [x for x in range(0, self.get_length()) if x not in rem]
        self.remove_columns(set(rem))
        #_LOG.debug("Alignment length reduced to %d" % len(subset))
        return subset

    def degap(self):
        """
        remove all
        """
        for name, seq in self.items():
            self[name] = self[name].replace("-", "")

    def get_hard_sub_alignment(self, sub_keys):
        """Creates a new alignment with a subset of the taxa."""
        new_alignment = MutableAlignment()
        new_alignment.datatype = self.datatype
        for key in sub_keys:
            new_alignment[key] = self[key]
        new_alignment.delete_all_gap()
        return new_alignment

class ReadonlySubalignment(ReadOnlyAlignment):
    """
    A class that emulates a subalignment of a given alignment. This is a
    readonly alignment and does not actually hold sequences in memory. It
    simply keeps a set of sequences that belong to the sub alignment,
    and implements methods of alignment class (and dictionaries) such that
    only subalignment keys are returned.
    """
    def __init__(self, keys, parent_alignment):
        self.seq_names = set(keys)
        self.parent_alignment = parent_alignment

    def __getitem__(self, key):
        if key in self.seq_names:
            return self.parent_alignment[key]
        else:
            raise KeyError(key)

    def __len__(self):
        return len(self.seq_names)

    def __iter__(self):
        for key in self.seq_names:
            yield key

    def get_mutable_alignment(self):
        ret = MutableAlignment()
        ret.set_alignment(self)
        return ret


class _AlignmentLookupHelper(object):
    """
    Internal helper calss
    """
    def __init__(self, pos, ref):
        self.pos = pos
        self.ref = ref

    def get(self, key, default=None):
        if key in self.ref:
            return self.ref[key][self.pos]
        else:
            return default

class ExtendedAlignment(MutableAlignment):
    """
    This is used to keep an extended alignment. An extended alignment
    has a bunch of sequences that are labeled as fragments. More importantly,
    columns of extended alignments are labeled with numbers and also it is
    known whether a column is an "insertion" column or a normal column.
    1) these can be read from .sto files.
    2) these alignments can be merged together.
    """
    def __init__(self, fragment_names):
        MutableAlignment.__init__(self)
        self.fragments = set(fragment_names)
        self._col_labels = []

    def set_alignment(self, alignment):
        MutableAlignment.set_alignment(self, alignment)
        self._reset_col_names()

    def read_file_object(self, file_obj, file_format='FASTA'):
        """ currently supports only fasta"""
        ret = MutableAlignment.read_file_object(self, file_obj, file_format)
        self._reset_col_names()
        return ret

    def remove_missing_fragments(self):
        for frag in list(self.get_fragment_names()):
            if frag not in self:
                self.fragments.remove(frag)

    def add_column(self, pos, char='-', new_label=None):
        """
        A new column is added to the alignment. The new label can be value,
        or one of three directives (see below).
        """
        MutableAlignment.add_column(self, pos, char)
        if new_label == "MAX":
            self._col_labels.insert(pos, max(self._col_labels) + 1)
        elif new_label == "INC_LAST":
            self._col_labels.append(max(self._col_labels) + 1)
        elif new_label == "RESET":
            self._reset_col_names()
        else:
            self._col_labels.insert(pos, new_label)

    def remove_column(self, pos, labels="REMOVE"):
        """
        Remove a column and potentially adjust column names.
        """
        MutableAlignment.remove_column(self, pos)
        if labels == "RESET":
            self._reset_col_names()
        elif labels == "REMOVE":
            self._col_labels = self._col_labels[:pos] + \
                self._col_labels[pos + 1:]

    def _get_col_labels(self):
        return self._col_labels
    col_labels = property(_get_col_labels)

    def _reset_col_names(self):
        """ sequentially label columns"""
        self._col_labels = list(range(0, self.get_length()))

    def get_fragments_readonly_alignment(self):
        """
        Return a readonly alignment that contains only the fragments.
        """
        return ReadonlySubalignment(self.fragments, self)

    def get_base_readonly_alignment(self):
        """
        Returns a readonly subalignment that does not contain fragments.
        """
        return ReadonlySubalignment(self.get_base_seq_names(), self)

    def get_fragment_names(self):
        return self.fragments

    def get_base_seq_names(self):
        return list(set(self.keys()) - self.fragments)

    def _read_sto(self, handle):
        """
        Reads a sto file and populates this current object. Figures out
        insertion columns by finding lower case letters, asterisks, and dots.
        """
        p = re.compile(r'[a-z*.]')
        insertions = []
        last_start_ind = -1
        for line in handle:
            line = line.strip()
            if line == "//":
                # End of the alignment. ignore meta-data
                break
            elif line == "" or line[0] == "#":
                pass
            else:  # not a comment
                key, seq = line.split()
                current = self.get(key, "")
                startind = len(current)
                # finding insertion columns in limited to first sequence
                if startind != last_start_ind:
                    # s = sum(len(x) for x in current)
                    insertions.extend([m.start() + startind
                                       for m in p.finditer(seq)])
                    last_start_ind = startind
                self[key] = current + seq.replace(".", "-")
#        for k,v in self.items():
#            self[k] = "".join(v)
        self._reset_col_names()
        return set(insertions)

    '''
    6.9.2022 - added by Chengze Shen
    a new function specifically dealing with sub-alignments of WITCH
    '''
    def read_query_alignment(self, query_name, path, aformat='fasta'):
        insertions = []
        if aformat == 'fasta':
            entries = [(n, s) for n, s in read_fasta(path)]
        elif isinstance(path, MutableAlignment):
            entries = [(n, s) for n, s in path.items()]

        num_elem_per_col = [0 for _i in range(len(entries[0][1]))]

        # count how many non-gaps in each col
        query_entry = None
        for i in range(0, len(entries), 1):
            name, entry = entries[i]
            if name != query_name:
                entry_count = tuple(1 if c != '-' else 0 for c in entry)
                num_elem_per_col = list(map(add, num_elem_per_col, entry_count))
            else:
                query_entry = entry
        if not query_entry:
            print('query entry {} missing in {}???'.format(query_name, path))

        # go over the query entry and see which column it has non-gap char
        # but backbone is gap (i.e., insertion)
        # additionally, record the columns that query aligns to
        query_entry = [c for c in query_entry]
        query_aligned_columns = []
        for j in range(0, len(query_entry), 1):
            if num_elem_per_col[j] == 0 and query_entry[j] != '-':
                query_entry[j] = query_entry[j].lower()
                insertions.append(j)
                query_aligned_columns.append(-1)
            elif num_elem_per_col[j] != 0 and query_entry[j] != '-':
                query_aligned_columns.append(j)
        self.fragments.add(query_name); self[query_name] = ''.join(query_entry)
        self._reset_col_names()

        # rename all columns to sequential numbers
        insertion = -1
        for c in insertions:
            self._col_labels[c] = insertion
            insertion -= 1
        regular = 0
        for c in range(0, self.get_length()):
            if self._col_labels[c] >= 0:
                self._col_labels[c] = regular
                regular += 1
        return query_name, set(insertions), tuple(query_aligned_columns)

    def build_extended_alignment(self, base_alignment, path_to_sto_extension,
                                 convert_to_string=True):
        """
        Given a base alignment and a path to an .sto file (or a list of paths),
        this methods populates self with an extended alignment by first reading
        the base alignment, and then merging in all the .sto extension
        alignments.
        Note that the .sto alignments should contain only fragments, and also
        note that there should be no taxon overlap among extensions, or between
        extensions and the base.
        """
        if isinstance(base_alignment, Alignment):
            self.set_alignment(deepcopy(base_alignment))
        elif isinstance(base_alignment, str):
            #_LOG.debug("Reading base alignment: %s." % base_alignment)
            self.read_filepath(base_alignment, "FASTA")

        if isinstance(path_to_sto_extension, str):
            paths = [path_to_sto_extension]
        else:
            paths = path_to_sto_extension

        self.from_string_to_bytearray()
        for path in paths:
            #_LOG.debug("Reading sto extended alignment: %s." % path)
            ext = ExtendedAlignment(self.fragments)
            ext.read_extended_alignment(path)
            #_LOG.info(
            #    "Merging extension sto file (%s) into base alignment (%s)." %
            #    (path, base_alignment))
            self.merge_in(ext, False)
            #_LOG.debug(
            #    ("Finished merging extension sto file (%s) into base "
            #     "alignment (%s).") %
            #    (path, base_alignment))
            del ext
        if convert_to_string:
            self.from_bytearray_to_string()

    def read_extended_alignment(self, path, aformat="stockholm",
                                assertion=False):
        """ Reads alignment from given path and figures out "insertion"
        columns. Labels insertion columns with special labels and labels the
        rest of columns (i.e. original columns) sequentially.
        """
        handle = open(path, 'r')
        insertions = None
        if aformat.lower() == "stockholm":
            insertions = self._read_sto(handle)
            #_LOG.debug("%s insertions: %d" % (path, len(insertions)))
        else:
            raise ValueError("format %s is not supported yet" % aformat)
        assert insertions is not None

        '''Assert that insertion columns have only gaps in original seqs and
        give them appropriate labels'''
        insertion = -1
        for c in insertions:
            if assertion:
                for k in self.get_base_seq_names():
                    assert not self[k][c] != "-", \
                        ("Insertion column has sequence among original "
                         "sequences. An error? column: %d k= %s" % (c, k))
            self.col_labels[c] = insertion
            insertion -= 1
        ''' Adjust labels of other columns '''
        i = 0
        for c in range(0, self.get_length()):
            if self.col_labels[c] >= 0:
                self.col_labels[c] = i
                i += 1

    def _is_insertion_label(self, i):
        """ This differentiates between an insertion and an original column
        """
        return i < 0

    def is_insertion_column(self, col):
        return self._is_insertion_label(self.col_labels[col])

    def relabel_original_columns(self, original_labels):
        """
        This methods relabels non-insertion columns in self based on the
        input labels. Insertion column labels will not be affected.
        """
        j = 0
        #_LOG.debug(
        #    "Relabeling %d (%d) with %d labels." %
        #    (self.get_length(), len(self._col_labels), len(original_labels)))
        for i in range(0, self.get_length()):
            if not self._is_insertion_label(self.col_labels[i]):
                assert \
                    j < len(original_labels), \
                    ("Not enough labels %d %d %d.\n %s" %
                        (i, j, len(original_labels), str(self._col_labels)))
                self.col_labels[i] = original_labels[j]
                j += 1
        assert\
            j == len(original_labels), \
            ("Some of original labels are unused."
             " Some columns from original alignment went missing? %d %d" %
                (j, len(original_labels)))

    def get_insertion_columns(self):
        return [i
                for (i, x)
                in enumerate(self.col_labels)
                if self._is_insertion_label(x)]

    def get_insertion_column_ranges(self):
        pos = []
        strt = None
        prev = None
        for p in self.get_insertion_columns():
            if prev is None or prev + 1 != p:
                if strt is not None:
                    pos.append((strt, prev))
                strt = p
            prev = p
        if strt is not None:
            pos.append((strt, prev))
        return pos

    def write_insertion_column_indexes(self, path):
        file_obj = open(path, 'w')
        r = ','.join("%d-%d" % (pair[0], pair[1])
                     for pair in self.get_insertion_column_ranges())
        file_obj.write(r)
        file_obj.write("\n")
        file_obj.close()

    def remove_insertion_columns(self):
        """
        Outputs a new alignment with insertion columns masked out.
        """
        cols = self.get_insertion_columns()
        s = []
        a = 0
        for b in cols:
            if b > a:
                s.append((a, b))
            a = b + 1
        s.append((a, len(self.col_labels)))
        for name, seq in list(self.items()):
            news = []
            for c in s:
                news.append(seq[c[0]:c[1]])
            self[name] = "".join(news)

    def write_insertion_maked_to_file(self, path):
        cols = self.get_insertion_columns()
        s = []
        a = 0
        for b in cols:
            if b > a:
                s.append((a, b))
            a = b + 1
        s.append((a, len(self.col_labels)))
        file_obj = open(path, 'w')
        for name, seq in self.items():
            file_obj.write('>%s\n' % name)
            for c in s:
                file_obj.write(seq[c[0]:c[1]])
            file_obj.write("\n")
        file_obj.close()

    def from_bytearray_to_string(self):
        for k, v in self.items():
            self[k] = v.decode()

    def from_string_to_bytearray(self):
        for k, v in self.items():
            self[k] = bytearray(v, encoding="utf8")

    def merge_in(self, other, convert_to_string=True):
        """
        Merges another alignment in with the current alignment. The other
        alignment needs to be an ExtendedAlignment as well. Since both
        alignments are extended alignments, we know column labels, and we also
        know which columns are insertions. Columns with the same labels are
        merged together. In other cases gaps are introduced as necessary to
        merge the two alignments. Insertion columns are considered different
        even when they have the same label.

        convert_to_string is by default true, meaning that in the beginning
        sequences of self are assumed to be string objects. These are converted
        to bytearray for better speed, but at the end of the merge, are
        converted back to string. When convert_to_string is False, no
        conversion back and from string is performed, *but* everything is
        assumed to be in bytearrays.
        So, self is assumed to be in bytearray before calling merge, and will
        remain in bytearray after calling merge. This is useful in cases where
        multiple alignments are merged in with self. Initially, everything
        should be turned into bytearray, and after all merging is finished,
        everything can be converted back to string (using
        from_bytearray_to_string and from_string_to_bytearray).
        """
        assert isinstance(other, ExtendedAlignment)
        #_LOG.debug("Merging started ...")
        if other.is_empty():
            return
        me = 0
        she = 0  # Assumption: alignments are female!
        me_len = self.get_length() if not self.is_empty() else 0
        she_len = other.get_length()
        insertion = -1

        merged_insertion_columns = 0

        ''' Add sequences from her to my alignment '''
        for f in other.fragments:
            self.fragments.add(f)
        if convert_to_string:
            self.from_string_to_bytearray()

        selfother = {}
        for k, v in other.items():
            # assert(k not in self,
            # "Merging overlapping alignments not implemented")
            if k not in self:
                selfother[k] = bytearray(v, encoding="utf8")
        while True:
            ''' Check exit conditions'''
            if me == me_len and she == she_len:
                break

            ''' Check the 5 possible statuses between she and I '''
            if she != she_len and other.is_insertion_column(she):
                if me != me_len and self.is_insertion_column(me):
                    ''' We both have a series of insertion columns'''
                    start = me
                    while(me != me_len and self.is_insertion_column(me) and
                          she != she_len and other.is_insertion_column(she)):
                        me += 1
                        she += 1
                        merged_insertion_columns += 1
                    run = me - start
                    self.col_labels[start:me] = list(range(
                        insertion, insertion-run, -1))
                else:
                    ''' Hers is a series of insertion columns'''
                    start = she
                    while she != she_len and other.is_insertion_column(she):
                        she += 1
                    run = she - start
                    ins = bytearray(b"-") * run
                    for seq in self.values():
                        seq[me:me] = ins
                    self._col_labels[me:me] = list(range(
                        insertion, insertion - run, -1))
                    insertion -= run
                    me += run
                    me_len += run
            elif me != me_len and self.is_insertion_column(me):
                ''' Mine is a series of insertion column'''
                start = me
                while me != me_len and self.is_insertion_column(me):
                    me += 1
                run = me - start
                ins = bytearray(b"-") * run
                for v in selfother.values():
                    v[start:start] = ins
                self.col_labels[start:me] = list(
                    range(insertion, insertion-run, -1))
                insertion -= run
            elif(she == she_len or (me != me_len and
                 self.col_labels[me] < other.col_labels[she])):
                ''' My column is not present (i.e. was allgap) in the
                    "other"'''
                start = me
                while(me < me_len and (she == she_len or me != me_len and
                      self.col_labels[me] < other.col_labels[she])):
                    me += 1
                run = me - start
                ins = bytearray(b"-") * run
                for v in selfother.values():
                    v[start:start] = ins
            elif(me == me_len or (she != she_len and
                 self.col_labels[me] > other.col_labels[she])):
                ''' Her column is not present (i.e. was allgap) in "me"'''
                start = she
                while(she < she_len and (me == me_len or she != she_len and
                      self.col_labels[me] > other.col_labels[she])):
                    she += 1
                run = she - start
                ins = bytearray(b"-") * run
                for seq in self.values():
                    seq[me:me] = ins
                self._col_labels[me:me] = other.col_labels[start:she]
                me += run
                me_len += run
            elif self.col_labels[me] == other.col_labels[she]:
                ''' A shared column'''
                while(me < me_len and she < she_len and
                      self.col_labels[me] == other.col_labels[she]):
                    she += 1
                    me += 1
            else:
                raise "hmmm, we thought this should be impossible? %d %d" % (
                    me, she)

        self.update(selfother)

        if convert_to_string:
            self.from_bytearray_to_string()
        #_LOG.debug("Merging finished ...")

        return merged_insertion_columns


################## Adapted from SEPP/jobs.py ####################
def readHMMSearch(path):
    assert isinstance(path, str), "Need to give me a path to HMMSearch result"
    assert os.path.isfile(path), "File path does not exist"
    f = open(path, 'r')
    lines = f.read().split('\n')

    results = {}
    pattern = re.compile(
    r"([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+"
    r"([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)")
    start_reading = False

    for line in lines:
        line = line.strip()
        # start reading when we see "E-value"
        if not start_reading and line.startswith("E-value") is True:
            start_reading = True
        # end reading if encounter empty line
        elif start_reading and line == "":
            start_reading = False
            break
        # extend the reading
        elif start_reading:
            matches = pattern.search(line)
            if matches is not None and matches.group(0).find("--") == -1:
                results[matches.group(9).strip()] = (
                        float(matches.group(1).strip()),
                        float(matches.group(2).strip()))
    f.close()
    return results

'''
Given an aligned string (represented by upper/lower case letters and gaps,
return a condensed version that have the lowercase letters from both sides
compressed to front/back.
'''
def compressInsertions(seq):
    p = re.compile(r'[A-Z]+')
    alns = [(m.start(), m.end()) for m in p.finditer(seq)]
    # do not perform such task if there is no aligned column at all
    if len(alns) == 0:
        return seq

    # first occurrence of aligned position defines the back of front search space
    # i.e., [start, end)
    f_start, f_end = 0, alns[0][0]
    f_len = f_end - f_start

    # last occurrence of aligned position defines the front of the back search space
    b_start, b_end = alns[-1][1], len(seq)
    b_len = b_end - b_start

    # simplest way of compression: remove all gaps and add them back
    f_str_ins = seq[f_start:f_end].replace('-', '')
    f_len_ins = len(f_str_ins)
    f_str = f_str_ins + '-' * (f_len - f_len_ins)

    b_str_ins = seq[b_start:b_end].replace('-', '')
    b_len_ins = len(b_str_ins)
    b_str = '-' * (b_len - b_len_ins) + b_str_ins

    # combine the compressed front/back with the remaining sequence
    combined = f_str + seq[f_end:b_start] + b_str

    return combined
