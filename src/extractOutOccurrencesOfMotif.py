import os 
import sys
from pbcore.io.FastaIO import FastaReader
import numpy

#motifs = sys.argv[1].split(",")
#modpositions = map(int, sys.argv[2].split(","))
genome = sys.argv[1]
median_position_fn = sys.argv[2]
motif_file=sys.argv[3]
tag=sys.argv[4]
motifs=[]
modpositions=[]

# for now recall motifs using a simple mapping to get rid of degenerate bases
# more generally we will need to actually account correctly for degenerate seq

def remap_bases (base):
	if base in "ACGT":
		return base
	else:
		return "N"

with open(motif_file) as f:
	for l in f:
		l = l.strip().split()
		print(l)
		m,i=l[0:2]
		print m
		m = "".join(map(remap_bases, list(m)))
		print m
		print i
		#m,i=l[0:1]
		motifs.append(m)
		modpositions.append(int(i))
		

posDict = {}
offsetDict = {}
count = 0
offset = 0
offsetlist = []
offsetcontig = []
#posToContig = {}
motifDict = {}
posToMotif = {}

def equalityCheck (seq1, seq2, m):
    # let seq1 be allowed to have degenerate bases
    for i in range(m):
        if not(seq1[i] == seq2[i] or seq1[i] == "N"):
            return False
    return True

for entry in FastaReader(genome):
    if count != 0:
        offset += 10
    seq = entry.sequence
    for i in range(len(motifs)):
        motif = motifs[i]
        modposition = modpositions[i]
        m  = len(motif)
        motifstr = motif+str(modposition)
        # check occurence of motif
        for i in range(len(seq)-len(motif)):
            # if dyad (i.e. ACTNNNNNNGAG) then do two equality checks:
            # 1.) seq[i:i+monad1] == dyad[0:monad1]
            # 2.) seq[i+monad1+lenN:i+monad1+lenN+monad2] == dyad[-monad2:]
            #if seq[i:i+m] == motif: # change for degenrate - make a function or a regex check
            if equalityCheck(motif, seq[i:i+m], m):
            #print entry.name, i+(modposition-1), motif
                pos = i+(modposition-1)
            #posDict.setdefault(entry.name, {})[pos]
                posDict[offset+pos] = 1
                posToContig = motifDict.setdefault(motifstr, {})
                posToContig[offset+pos] = entry.name
                posToMotif[offset+pos] = motifstr
    offsetDict[entry.name] = offset
    offset += len(entry.sequence)
    offsetlist.append(offset)
    offsetcontig.append(entry.name)
    count += 1

# open median position and extract ipds
#contigIpds = {}
motifIpds = {}
with open (median_position_fn) as f:
    for l in f:
        ll = l.strip().split()
        if len(ll) < 3:
            break
        pos, ipd, strand = ll
        if strand == "0": # quick hack we're ignoring reverse strand
            continue
        pos = int(pos)
        ipd = float(ipd)
        if pos in posDict:
            motif = posToMotif[pos]
            contigIpds = motifIpds.setdefault(motif, {})
            posToContig = motifDict[motif]
	    pos_ipd = [pos,ipd]
	    tmp = str(pos_ipd)[1:-1]
	    print tmp 
            contigIpds.setdefault(posToContig[pos], []).append(tmp)


# contig, mean IPD, motif 
on="%s.cim" % tag
out_name=open(on, 'w')
# TBD print out full dataset in a table
motifkeys = motifIpds.keys()
motifkeys.sort()
for motif in motifkeys:
    contigIpds = motifIpds[motif]
    for contig, values in contigIpds.items():
	#print >>sys.stderr, motif, contig, numpy.median(values), numpy.mean(values)
        for val in values:
		m=motif[:-1]
		#print tmp_val
		out_name.write("%s %s %s\n" % (contig, val, m))



