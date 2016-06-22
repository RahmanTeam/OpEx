# Sort BAM file

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + '/../../pysamdir')
import pysam

pysam.sort("-",sys.argv[1])
