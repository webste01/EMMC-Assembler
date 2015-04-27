#get the ratio between either the mean or median native data and predicted ipds
import os
import sys
from collections import defaultdict
import re
import numpy as np

name=sys.argv[1]
nat=open(sys.argv[1],'r')
pred=open(sys.argv[2],'r')
fn_ratio="%s.ratio" % name
f_ratio = open(fn_ratio,"w")
#fs=defaultdict(list)
#ns=defaultdict(list)
fs={}
ns={}
print >>sys.stderr, "reading native"
linecount = 0
for line in nat:
	line=line.split()
	pos=int(line[0])
	ipd=float(line[1])
	strand=line[2]
	if strand == "1":
		fs[pos]=ipd
	else:
		ns[pos]=ipd
	if linecount % 1000000 == 0:
		print >>sys.stderr, linecount
	linecount += 1
#fs_pred=defaultdict(list)
#ns_pred=defaultdict(list)

fs_pred={}
ns_pred={}
print >>sys.stderr, "reading in silico"
for line in pred:
	line=line.split()
	pos=int(line[0])
	ipd=float(line[1])
	strand=line[2]
	if strand == "1":
		fs_pred[pos]=ipd
	else:
		ns_pred[pos]=ipd
	if linecount % 1000000 == 0:
		print >>sys.stderr, linecount
	linecount += 1

print >>sys.stderr, "printing ratios"
for pos in fs:
	if pos % 1000000 == 0:
		print >>sys.stderr, pos

	if pos in fs_pred:
		print pos
		print fs[pos]
		print fs_pred[pos]	
		rat=np.divide(fs_pred[pos],fs[pos])
		print rat
		strand=1
		f_ratio.write("%s %.5f %s\n" % (pos, rat, strand))
for pos in ns:
	if pos % 1000000 == 0:
		print >>sys.stderr, pos

        if pos in ns_pred:
		print pos
		print ns[pos]
		print ns_pred[pos]
        	rat=np.divide(ns_pred[pos],ns[pos])
        	print rat
		strand=0
        	f_ratio.write("%s %.5f %s\n" % (pos, rat, strand))
f_ratio.close()
