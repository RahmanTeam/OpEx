# Converting SAM file to uncompressed BAM file

import pysam
infile = pysam.AlignmentFile("-", "r")
outfile = pysam.AlignmentFile("-", "wbu", template=infile)
for s in infile: outfile.write(s)

