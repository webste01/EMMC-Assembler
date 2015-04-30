#!/bin/bash
############
# BASH script mastercode:
############
REF=$1
CMPH5=$2
motifs=$3
best_edges=$4
tag=$5
gkpStore=$6
tigStore=$7
src=src

#convert celera assembler best edges file to gml
python ${src}/CA_best_edge_to_GML_ali.py ${gkpStore} ${tigStore} ${best_edges} $
{best_edges}.gml

#Get in silico IPDs. Output is ${REF}.pred_ipds
python ${src}/get_in_silico_ipds.py ${REF} 

#CMPH5 normalization 
#outputs nativedata.txt
python ${src}/cmph5SubreadNormMultiRef.py ${REF} ${CMPH5} 

#Get reduced IPDs mean and median **bake into one script
#outputs nativedata.txt.reduced_mean nativedata.txt.reduced_median
python ${src}/med_mean_ipd_calc.py nativedata.txt

#get insilico/native reduced ratio
python ${src}/get_ratio.py native.txt.reduced_mean ${REF}.pred_ipds
#python ${src}/get_ratio.py native.txt.reduced_median ${REF}.pred_ipds

#extract motifs
python ${src}/extractOutOccurrencesOfMotif.py ${REF} nativedata.txt.reduced_mean
 ${motifs} ${tag}
#get gml from celera assembler
python ${src}/CA_best_edge_to_GML_ali.py ${best_edges} 

#graphing and statistical testing
python ${src}/graph_coloring_cut_ipds.py ${best_edges}.gml ${tag}.cim ${tag} ${R
EF}
