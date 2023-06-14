# EMMA - Extending Multiple alignments with MAFFT--add
<a href="https://archive.softwareheritage.org/browse/origin/?origin_url=https://github.com/c5shen/EMMA">
    <img src="https://archive.softwareheritage.org/badge/origin/https://github.com/c5shen/EMMA/" alt="Archived | https://github.com/c5shen/EMMA"/>
</a>


(C) Chengze Shen

EMMA is an ensemble usage of `MAFFT --add` (particularly, `MAFFT` with `-linsi` option) on large datasets. On the MAFFT webpage, `MAFFT-linsi --add` is accurate for adding sequences to an existing alignment, but is only recommended to use for a few hundreds of sequences. This project aims to scale `MAFFT-linsi --add` to run on large datasets with several thousands of sequences with similar alignment accuracy.

----
News
----
1. Currently developing an extension to the accepted version of EMMA at WABI 2023. Now support the assignment and alignment of query sequences to happen at different subsets. That is, we use smaller/less diverse subsets to assign query sequences (could be more accurate for assignment), but use the corresponding larger subsets (the subset that decomposes to the assigned smaller subsets) to align the query sequences (presumably more accurate due to inclusion of more query sequences for `MAFFT-linsi--add`).


---------------
Method Overview
---------------
### Algorithm
Given an input existing alignment $C$ on set $S$ (i.e., backbone alignment) and a set of unaligned sequences $Q$ (i.e., query sequences) that we want to add to $C$, EMMA outputs an alignment on $S\cup Q$ that induces $C$ when restricted to $S$. The detailed pipeline is presented below:
1. __Construct a set of constraint sub-alignments from $C$__: Decompose $C$ to sub-alignments using the [UPP](https://github.com/smirarab/sepp/blob/master/README.UPP.md) decomposition strategy but limit to sub-alignments with $|Q_i|$ sequences, $50\leq |Q_i|\leq u$ ($u$ a user-provided free parameter).
2. __Define the set of sub-problems__: Assign each query sequence $q\in Q$ to the best-fitting sub-alignment from Step 1. The assignment is determined by first constructing HMMs on the sub-alignments and then selecting the HMM with the highest bitscore for each $q$.
3. __Run `MAFFT-linsi--add` on each sub-problem__: For each sub-problem (i.e., a sub-alignment $C_i$ on set $S_i$ with assigned query sequences $Q_i$), construct an extended sub-alignment on $S_i\cup Q_i$ using `MAFFT-linsi--add`.
4. __Merge extended sub-alignments using transitivity__: All extended sub-alignments are consistent with each other (see proof in the main paper) and can merge to the backbone alignment with transitivity (see [SEPP/UPP](https://github.com/smirarab/sepp)). The merging produces the final alignment on $S\cup Q$.

#### Publication
* Currently accepted in WABI 2023.

------------
Installation
------------
EMMA was tested and benchmarked on the following systems:
* Red Hat Enterprise Linux Server release 7.9 (Maipo) with __Python 3.7.0__
* Ubuntu 22.04 LTS with __Python 3.7.12__

EMMA requires the usage of `MAFFT` binaries. One is provided with the package (v7.490 2021/Oct/30), but the `MAFFT` binaries in the user's `$PATH` environment are prioritized. If you experience any difficulties running EMMA, please contact Chengze Shen (chengze5@illinois.edu).

### Requirements
```
python>=3.7
configparser>=5.0.0
dendropy>=4.5.3
numpy>=1.15
psutil>=5.0
scipy>=1.1.0
```

### Installation Steps
```bash
# 1. Install via GitHub repo
git clone https://github.com/c5shen/EMMA.git

# 2. Install all requirements
cd EMMA
pip3 install -r requirements.txt

# 3. Use emma.py, -h to see allowed commandline parameters
python3 emma.py [-h]
```

-------
Examples
-------
### Scenario A: given an input alignment, adding a set of unaligned sequences
```bash
python3 emma.py -b [input alignment] -q [unaligned sequences] -d [output directory] -o est.aln.fasta
```

# TO-DO
* Add more customizable configuration support as WITCH.
* ~Finish up the pipeline so it supports building an alignment from scratch and not relying on UPP output.~
