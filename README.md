# EMMC-Assembler
EMMC-Assembler
========

**E**pigenetic **M**etagenomic **M**otif **C**haracterization **Assembler**

## Overview
This tool compares the distribution of IPDs between NATIVE and WGA data for every existing kmer in the genome (using log-likelihood rations) and then uses an iterative greedy procedure to identify significant motifs.  The output of this tool is a list of the significantly modified kmers.

##Running
To run the code simply type:

<code>sh EMMCAssemble.sh [reference.fasta] [nat.cmp.h5] [motif file] [best.edges] [output directory] </code>

### INPUTS

reference.fasta: reference fasta file for the metagenomic sample
nat.cmp.h5: cmph5 file for the native (NOT WGA) metagenomic sample
motif file: a file listing the motifs and epigenetically modified positions. An example is included here as "example_motif.txt"
best.edges: the best.edges file produced by Celera Assembler in the tigstore

### INSTALLATION/REQUIREMENTS

<ul>
<li>PYTHON </li>
<li>pbtools </li>
<li>pbcore</li>
<li></li>
<li></li>
<li></li>
<li></li>
<li></li>
<li></li>
</ul>
 
## Citation
To cite the code, please refer to:

##Questions
Please contact elizabeth.webster@icahn.mssm.edu
