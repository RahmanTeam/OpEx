# Sort BAM file

import sys
import pysam
pysam.sort("-",sys.argv[1])