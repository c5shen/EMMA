# EMMA - Extending Multiple alignments with MAFFT-linsi --add
<a href="https://archive.softwareheritage.org/browse/origin/?origin_url=https://github.com/c5shen/EMMA">
    <img src="https://archive.softwareheritage.org/badge/origin/https://github.com/c5shen/EMMA/" alt="Archived | https://github.com/c5shen/EMMA"/>
</a>


(C) Chengze Shen

EMMA is an ensemble usage of `MAFFT --add` (particularly, `MAFFT` with the `-linsi` option) on large datasets. On the MAFFT webpage, `MAFFT-linsi --add` is accurate for adding sequences to an existing alignment but is only recommended for a few hundred of sequences. This project aims to scale `MAFFT-linsi --add` to run on large datasets with hundreds of thousands of sequences with similar (and sometimes better) alignment accuracy.

----
News
----
2. (NEW) Now will not exit after creating `main.config` for the first time (if directly running `emma.py [arguments]` without running `setup.py`). In addition, now EMMA will detect if the binaries defined in `main.config` are runnable and notify users if some binaries have an error when executing.
1. Checkpoint system! Now you can resume from any point if a previous run was interrupted somehow (except for the HMMSearch step, currently in implementation).
2. Now automatically detects input data type/molecule (`amino`, `dna`, or `rna`).
3. Now has a progress bar for all intermediate steps (for better progress tracking!).


# TO-DO
* ~Add checkpoint support.~ Still need HMMSearch step checkpoint system.
* ~Add more customizable configuration support as WITCH.~
* ~Finish up the pipeline so it supports building an alignment from scratch and not relying on UPP output.~


---------------
Method Overview
---------------
### Algorithm
Given an input existing alignment $C$ on set $S$ (i.e., backbone alignment) and a set of unaligned sequences $Q$ (i.e., query sequences) that we want to add to $C$, EMMA outputs an alignment on $S\cup Q$ that induces $C$ when restricted to $S$. The detailed pipeline is presented below:
1. __Construct a set of constraint sub-alignments from $C$__: Decompose $C$ to sub-alignments using the [UPP](https://github.com/smirarab/sepp/blob/master/README.UPP.md) decomposition strategy but limit to sub-alignments with $|Q_i|$ sequences, $l\leq |Q_i|\leq u$ ($l,u$ are user-provided free parameters; default values are $l=10,u=25$). This step creates a set of subsets that can overlap.
2. __Define the set of sub-problems__: Assign each query sequence $q\in Q$ to the best-fitting sub-alignment from Step 1. The assignment is determined by first constructing HMMs on the sub-alignments and then selecting the HMM with the highest adjusted bitscore (see [WITCH](https://github.com/c5shen/WITCH)) for each $q$.
3. __Run `MAFFT-linsi--add` on each sub-problem__: For each sub-problem (i.e., a sub-alignment $C_i$ on set $S_i$ with assigned query sequences $Q_i$), construct an extended sub-alignment on $S_i\cup Q_i$ using `MAFFT-linsi--add`.
4. __Merge extended sub-alignments using transitivity__: All extended sub-alignments are consistent with each other (see proof in the main paper) and can merge to the backbone alignment with transitivity (see [SEPP/UPP](https://github.com/smirarab/sepp)). The merging produces the final alignment on $S\cup Q$.

#### Publication
* Currently accepted in WABI 2023.
* Currently published on Algorithms of Molecular Biology (https://doi.org/10.1186/s13015-023-00247-x).

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
dendropy>=4.5.2,<4.6.0
numpy>=1.15
psutil>=5.0
scipy>=1.1.0
tqdm>=4.0.0
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

#### IMPORTANT: macOS users
I currently do not have a working macOS system at hand for software testing purposes, hence the compiled macOS binaries included in EMMA may not be executable right now (reported by some users).
    
* The most likely issue right now should be `FastTreeMP` (`main.config`, defined under `[Basic]` and `[MAGUS]` as `fasttreepath`). In the case this binary file is not working out, please download the source code from [FastTree.c](http://www.microbesonline.org/fasttree/FastTree.c).
Direct quote from the [FastTree 2.1 webpage](http://www.microbesonline.org/fasttree/#Install) below.
* After compilation, you can either put `FastTreeMP` under your `$PATH` variable, or change `fasttreepath` in `main.config` to point to your compiled FastTree executable.
* The other possible culprit `mafft`, for which you can download or compile the executable from [MAFFT for macOS official page](https://mafft.cbrc.jp/alignment/software/macosx.html). You can follow the same step above to change the `mafftpath` in `main.config`.

Quote for FastTree:

> If you use a Mac or other platform not included above, download FastTree.c and run
>
> `gcc -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm`
>
> (`gcc` is installed on many Mac OS X and Unix machines. If you use a Mac, you may need to install it from Xcode. `gcc` is also available for virtually every platform.) Note that FastTree will try to use SSE2/SSE3 instructions to speed up some inner loops. This will not work on many Windows or Mac machines. If FastTree will not run, then try compiling it with this command instead:
>
> `gcc -DNO_SSE -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm`
>
> We have also heard that the -finline-functions option can cause an error. You can omit this option.
>
> If you want to build the multi-threaded "FastTreeMP," use
>
> `gcc -DOPENMP -fopenmp -O3 -finline-functions -funroll-loops -Wall -o FastTreeMP FastTree.c -lm`

-------
Examples
-------
Scripts of the following examples can be found in `example/run.sh`. You can run each scenario with
```bash
./run.sh [i]    # i can be 1, 2, or 3
```

### Scenario 1: given an input alignment and its tree, add unaligned sequences
```bash
python3 emma.py -b [input alignment] -e [input tree] -q [unaligned sequences] -d [output directory] -o est.aln.fasta
```

### Scenario 2: given just an input alignment, add unaligned sequences
```bash
python3 emma.py -b [input alignment] -q [unaligned sequences] -d [output directory] -o est.aln.fasta
```

### Scenario 3: given just unaligned sequences, align them all
```bash
# > the "backbone sequences" will be selected from inputs and aligned with default MAGUS
# > a tree will be created for the backbone alignment using FastTree2
python3 emma.py -i [input sequences] -d [output directory] -o est.aln.fasta
```
