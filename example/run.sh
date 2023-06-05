#!/bin/bash
#SBATCH --time=03:59:59
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --partition=cs
#SBATCH --mem=32GB
#SBATCH -o %j.out
#SBATCH -e %j.err

module load python
module load java/11
time=/usr/bin/time

emafftaddbin=../emma.py
t=8

query_path=./data/queries.unaln.fasta
backbone_path=./data/backbone.aln.fasta
backbone_tree=./data/backbone.est.tre
ref_path=./data/queries.aln.fasta

molecule=dna
lower=10
upper=30
alignment_size=200
subproblem_size=500

outdir=./example_output
if [ ! -d $outdir ]; then
    mkdir -p $outdir
fi

if [ ! -s $outdir/est.aln.fasta ]; then
    rm -r $outdir; mkdir -p $outdir
    { $time -v python3 $emafftaddbin -t $t -b ${backbone_path} \
        -e ${backbone_tree} -q ${query_path} -d $outdir \
        --molecule $molecule \
        --lower $lower \
        --upper $upper \
        --alignment-size ${alignment_size} \
        --subproblem-size ${subproblem_size} \
        -o est.aln.fasta ; } 2> $outdir/runtime.txt
fi
