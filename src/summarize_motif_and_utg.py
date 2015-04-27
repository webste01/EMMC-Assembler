import os
import sys
import numpy as np
cimfile = sys.argv[1]

motifUtgDict = {}
utgMotifDict = {}
with open (cimfile) as f:
    for l in f:
        ll = l.strip().split()
        utg, val, mot = ll
        val = float(val)
        motifUtgDict.setdefault(mot, {}).setdefault(utg, []).append(val)
        utgMotifDict.setdefault(utg, {}).setdefault(mot, []).append(val)

'''
for mot, utgDict in motifUtgDict.items():
    for utg, vals in utgDict.items():
        print utg, np.mean(vals), mot
        #print utg, 
        #print utg,
'''
motifs = motifUtgDict.keys()
print "utg motif1 motif2 motif3 motif4" 
for utg, motDict in utgMotifDict.items():
    print utg,
    for mot in motifs:
        if mot in motDict:
            vals  = motDict[mot]
            print np.mean(vals),
        else:
            print 0,
    print
            

