#cd /sc/orga/scratch/webste01/in_silico/kineticsTools/src 
import pbtools
import pbcore
import logging
import os
import sys
import re
import h5py
import numpy as np
import ctypes as C
from pbtools.kineticsTools.sharedArray import SharedArray
from pbtools.kineticsTools.ipdModel import IpdModel
from pbcore.io.FastaIO import FastaReader
ref_size=[]
ref=sys.argv[1]
ref_len_dict = {}
offsetDict = {}
offset = 0
for entry in FastaReader(ref):
	a = len(entry.sequence)
	ref_size.append(a)
	ref_len_dict[entry.name] = a
	offsetDict[entry.name] = offset
	offset += a
#	ref_size=len(entry.sequence)
# Only takes the first contig! - now takes in all the contigs (4)
#break
print ref_size
print len(ref_size)
"""
Generate a barebones reference entry that this kinetics code needs:
(1) contig.sequence
(2) contig.id
"""
def create_fastaRecord( fasta ):
    class Contig:
            def __init__( self, entry ):
                    self.id=entry.name
                    self.sequence=entry.sequence
    contigs_diff=[]
    for entry in FastaReader(fasta):
            contig=Contig(entry)
            contigs_diff.append(contig)
    return contigs_diff
"""
Use the surrounding sequence context of the reference positions to predict
IPD values to be used as the control.


# FIXME!!!!!!!!
******* This will initially only support single contig references. *******
"""
contigs= create_fastaRecord( ref )
#lut file from: https://github.com/khaikhai/kineticsTools/blob/master/src/pbtools/kineticsTools/kineticLut.h5
LUTh5 = "/sc/orga/scratch/webste01/in_silico/kineticsTools/src/pbtools/kineticsTools/kineticLut.h5"
model = IpdModel( fastaRecords=contigs, modelFile=LUTh5 )
################
# This is where it needs some updates/love. Specifically lines 395-410 in the pbtools.kineticsTools.ipdModel
# file.
################


fn= "%s.pred_ipds" % ref
f= open(fn, "w")
ref_pred_IPDs = {0:{}, 1:{}}
#Will need to adjust for multiple contigs
for strand in ref_pred_IPDs.keys():
	for ref_contig_id, contig_length in ref_len_dict.items():
		#ref_contig_id   = model.refDict.keys()[0]
		meanIpdFunc     = model.predictIpdFunc( refId=ref_contig_id )
		cognateBaseFunc = model.cognateBaseFunc(ref_contig_id)
		print ref_contig_id
		print ref_size.index(contig_length)
		print contig_length 
		print strand
		for pos in range(contig_length): 
			ref_base = cognateBaseFunc(pos,strand)
			ipd= meanIpdFunc( pos, strand )
			# Make sure we stay in the 0-based position space!
			ref_pred_IPDs[strand][pos] = ipd
			newpos = pos + offsetDict[ref_contig_id]
   			f.write("%i %.5f %i\n" % (newpos, ipd, strand))
f.close()

