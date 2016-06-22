# ...

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + '/../../pysamdir')
import pysam

pysam.index(sys.argv[1]+'.bam')
