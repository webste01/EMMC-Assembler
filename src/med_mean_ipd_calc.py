#Get the mean IPD value
import os
import sys
from collections import defaultdict
import re
import numpy as np
nat=open(sys.argv[1],'r')
name=sys.argv[1]
fs=defaultdict(list)
ns=defaultdict(list)
for line in nat:
	line=line.split()
	pos=line[0]
	ipd=float(line[1])
	strand=line[2]
	if strand == "1":
		fs[pos].append(float(ipd))
	else:
		ns[pos].append(float(ipd))
fn_med="%s.reduced_median" % name
fn_mean="%s.reduced_mean" % name
f_med= open(fn_med,"w")
f_mean=open(fn_mean,"w")
median_fs=defaultdict(list)
median_ns=defaultdict(list)
mean_fs=defaultdict(list)
mean_ns=defaultdict(list)

for pos in fs.keys():
	tmp = [ float(ipd) for ipd in fs[pos] ]	
	median_tmp=np.median(tmp)
	median_fs[pos]=median_tmp
	mean_tmp=np.mean(tmp)
        mean_fs[pos]=mean_tmp
	strand=1
	f_med.write("%s %.5f %s\n" % (pos, median_tmp, strand))
        f_mean.write("%s %.5f %s\n" % (pos, mean_tmp, strand))
for pos in ns.keys():
	tmp = [ float(ipd) for ipd in ns[pos] ]
	median_tmp=np.median(tmp)
	median_ns[pos]=median_tmp
        mean_tmp=np.mean(tmp)
        mean_ns[pos]=mean_tmp
	strand=0
        f_med.write("%s %.5f %s\n" % (pos, median_tmp, strand))
        f_mean.write("%s %.5f %s\n" % (pos, mean_tmp, strand))
f_med.close()
f_mean.close()

