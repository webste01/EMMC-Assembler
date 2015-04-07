# EMMC-Assembler

========

**E**pigenetic **M**etagenomic **M**otif **C**haracterization **Assembler**

## Overview
This tool uses epigenetically modified motifs to fine map the bacterial genomes present in mixed metagenomic samples. The output is fully assembled genomes present in the metagenomic sample. 

##Running
To run the code simply type:

<code>sh EMMCAssemble.sh [reference.fasta] [nat.cmp.h5] [motif file] [best.edges] [name] </code>

### INPUTS

reference.fasta: reference fasta file for the metagenomic sample

nat.cmp.h5: cmph5 file for the native (NOT WGA) metagenomic sample

motif file: a file listing the motifs and epigenetically modified positions. An example is included here as "example_motif.txt"

best.edges: the best.edges file produced by Celera Assembler from the unitigger

name: the name tag for the files created

### INSTALLATION/REQUIREMENTS

<ul>
<li>PYTHON </li>
<li>pbtools </li>
<li>pbcore</li>
<li>numpy</li>
<li>h5py</li>
<li>shlex</li>
<li>subprocess</li>
<li>networkx</li>
<li>matplotlib</li>
</ul>
 
## Citation
To cite the code, please refer to:

##Questions
Please contact elizabeth.webster@icahn.mssm.edu
