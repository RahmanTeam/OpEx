# Converting SAM file to uncompressed BAM file

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + '/../../pysamdir')
import pysam

infile = pysam.Samfile("-", "r")
outfile = pysam.Samfile("-", "wbu", template=infile)

for s in infile: outfile.write(s)


