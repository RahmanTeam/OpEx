#!/bin/bash

# Retrieving directory and Python version
dir=$(pwd)
pver=$(python -c "import sys; print sys.version[:3]")

# Creating temporary and final directories for pysam
mkdir -p $dir/tmp/lib/python$pver/site-packages/
mkdir $dir/pysamdir

# Unpacking pysam 0.7.7
tar -zxvf pysam-0.7.7.tar.gz

# Building pysam 0.7.7
cd pysam-0.7.7
python setup.py build

# Installing pysam 0.7.7 to temporary directory
PYTHONPATH="${PYTHONPATH}:$dir/tmp/lib/python$pver/site-packages"
export PYTHONPATH 
python setup.py install --prefix=$dir/tmp
cd ..

# Copying pysam from temporary to final directory
dest=$(ls -d $dir/tmp/lib/python2.7/site-packages/*/)
cp -r $dest/* $dir/pysamdir/

# Removing temporary directories of pysam
rm -r tmp
rm -r pysam-0.7.7


