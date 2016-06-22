import os
import sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + '/../../pysamdir')
import pysam

print reduce(lambda x, y: x + y, [ eval('+'.join(l.rstrip('\n').split('\t')[2:]) ) for l in pysam.idxstats(sys.argv[1]) ])
