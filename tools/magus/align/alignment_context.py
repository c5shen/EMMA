'''
Created on Dec 4, 2020

@author: Vlad
'''

import os
from helpers import sequenceutils
from tasks import task, manager
from configuration import Configs

'''
The AlignmentContext data structure maintains all information pertaining to a single MAGUS alignment.
The main thread keeps one context active at a time (to manage resources). Subalignments will spawn their
own contexts, which will become active when the parent context is in a blocking wait. The previously active
context resumes when the current alignment is completed.
'''

class AlignmentContext:
        
    def __init__(self, **kwargs):
        self.outputFile = None
        self.workingDir = None
        self.sequencesPath = None
        self.subsetPaths = []
        self.subalignmentPaths = []
        self.backbonePaths = []
        self.guideTree = None
        
        self.unalignedSequences = None
        #self.taxa = []
        self.subsets = []
        self.subalignments = []
        self.taxonSubsetMap = {}
        self.taxonSubalignmentMap = {}
                
        self.backboneTaxa = {}
        self.backboneExtend = set()
        self.backboneSubalignment = {}
        
        self.subalignmentTasks = []
        self.backboneTasks = []
        self.graph = None
        
        '''
            6.1.2021 - added by Chengze Shen
            a dict: <key> = backbone aln name, <value> = weights
            a path: the file path to define the dict
        '''
        self.backboneWeights = {}
        self.backboneWeightsPath = None

        for attr in kwargs:
            vars(self)[attr] = kwargs.get(attr)
        
        if not os.path.exists(self.workingDir):
            os.makedirs(self.workingDir)
    
    def awaitSubalignments(self):
        task.awaitTasks(self.subalignmentTasks)
    
    def initializeSequences(self):
        self.unalignedSequences = {}
        for i, subsetPath in enumerate(self.subsetPaths):
            self.subsets.append([])
            subset = sequenceutils.readFromFastaOrdered(subsetPath, removeDashes=True)
            for sequence in subset:
                self.unalignedSequences[sequence.tag] = sequence
                self.taxonSubsetMap[sequence.tag] = i
                self.subsets[i].append(sequence.tag)
        
        if Configs.constrain:
            self.subalignments = self.subsets
            self.taxonSubalignmentMap = self.taxonSubsetMap
        else:
            for s in self.subsets:
                for taxon in s:
                    self.taxonSubalignmentMap[taxon] = len(self.subalignments)
                    self.subalignments.append([taxon])
       
    '''
        6.1.2021 - added by Chengze Shen
        a new function to initialize readings of weights from the
        backbone weights path (if such path exists)
    '''
    def initializeBackboneWeights(self):
        if self.backboneWeightsPath:
            # the weights should be put in the following manner:
            # > each line denotes a weight (for a backbone)
            # > for each line, it should have the format: backbone path,weight
            with open(self.backboneWeightsPath, 'r') as f:
                lines = f.read().split('\n')[:-1]
                for line in lines:
                    if line == '':
                        continue
                    w = [x.strip() for x in line.split(',')]
                    assert len(w) == 2
                    self.backboneWeights[w[0]] = float(w[1])
        else:
            return


    def initializeBackboneSequenceMapping(self):
        if len(self.backboneTaxa) == 0:
            backboneSubsetTaxonMap = {i : subset for i, subset in enumerate(self.subsets)}
        else:
            backboneSubsetTaxonMap = {}
            for taxon in self.backboneTaxa:
                i = self.taxonSubsetMap[taxon]
                backboneSubsetTaxonMap[i] = backboneSubsetTaxonMap.get(i, [])
                backboneSubsetTaxonMap[i].append(taxon) 
        
        if Configs.constrain:
            for i, subalignPath in enumerate(self.subalignmentPaths):
                subalignment = sequenceutils.readFromFasta(subalignPath, removeDashes=False)
                for taxon in backboneSubsetTaxonMap.get(i, []):
                    self.backboneSubalignment[taxon] = subalignment[taxon]
        else:
            self.backboneSubalignment = self.unalignedSequences
    
    def __enter__(self):
        manager.TaskManager.contextStack.append(self)
        return self
            
    def __exit__(self, excType, excVal, excTb):
        manager.TaskManager.contextStack.pop()
    
