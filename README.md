# eMAFFTadd
eMAFFTadd is an ensemble usage of `MAFFT-linsi --add` on large datasets. On the MAFFT webpage, `MAFFT-linsi --add` is only recommended to a few hundreds of sequences. This project aims to scale `MAFFT-linsi --add` to run on large datasets with several thousands of sequences without losing accuracy.

_Currently the software requires UPP HMMsearch output files as inputs, since it does not calculate bit-scores between sequences and HMMs itself._

# TO-DO
* finish up the pipeline so it supports building an alignment from scratch and not relying on UPP output.

# Usage

# Example
