#!/bin/bash
#SBATCH --time=03:59:59
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --partition=cs
#SBATCH --mem=32GB
#SBATCH -o %j.out
#SBATCH -e %j.err

#module load python
#module load java/11
time=/usr/bin/time

emafftaddbin=../emma.py
t=-1

scenario=1
if [[ $1 != "" ]]; then
    scenario=$1
fi

query_path=./data/queries.unaln.fasta
backbone_path=./data/backbone.aln.fasta
backbone_tree=./data/backbone.est.tre
ref_path=./data/queries.aln.fasta
unaln_path=./data/all.unaln.fasta

lower=10
upper=25
alignment_size=$lower
subproblem_size=500

# scenario 1 - both backbone alignment and tree are available
if [[ $scenario == 1 ]]; then
    outdir=scenario_1_output
    python3 $emafftaddbin -t $t -b ${backbone_path} -e ${backbone_tree} \
        -q ${query_path} -d $outdir -o est.aln.fasta \
        --keep-decomposition
# scenario 2 - backbone alignment available but backbone tree is missing
elif [[ $scenario == 2 ]]; then
    outdir=scenario_2_output
    python3 $emafftaddbin -t $t -b ${backbone_path} -q ${query_path} \
        -d $outdir -o est.aln.fasta \
        --keep-decomposition
# scenario 3 - both backbone alignment and tree are missing; only input sequences
elif [[ $scenario == 3 ]]; then
    outdir=scenario_3_output
    python3 $emafftaddbin -t $t -i ${unaln_path} -d $outdir -o est.aln.fasta \
        --keep-decomposition
# scenario 4 - use experimental setting for query assignment and alignment
elif [[ $scenario == 4 ]]; then
    outdir=scenario_4_output
    python3 $emafftaddbin -t $t -b ${backbone_path} -e ${backbone_tree} \
        -q ${query_path} -d $outdir -o est.aln.fasta \
        --experimental --lower 10 --upper 25 --alignment-size 25 \
        --keep-decomposition
fi
