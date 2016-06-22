#!/bin/bash

message () {
echo; printf '%.0s=' {1..80}; echo
echo "OpEx: $1"
printf '%.0s=' {1..80}; echo; echo
}

# Current directory
opexdir=$(pwd)

# Download BWA
message "Downloading BWA" 
cd tools
wget http://sourceforge.net/projects/bio-bwa/files/bwa-0.5.10.tar.bz2 --no-check-certificate
tar -jxvf bwa-0.5.10.tar.bz2
rm bwa-0.5.10.tar.bz2; echo

# Download Stampy
message "Downloading Stampy"
wget http://www.well.ox.ac.uk/bioinformatics/Software/stampy-v1.0.14.1.tgz
tar zxvf stampy-v1.0.14.1.tgz
rm stampy-v1.0.14.1.tgz; echo

# Download Picard
message "Downloading Picard"
wget http://sourceforge.net/projects/picard/files/picard-tools/1.48/picard-tools-1.48.zip --no-check-certificate
unzip picard-tools-1.48.zip
rm picard-tools-1.48.zip; echo

# Download Platypus
message "Downloading Platypus"
wget http://www.well.ox.ac.uk/bioinformatics/Software/Platypus_0.1.5.tgz
tar zxvf Platypus_0.1.5.tgz
mv PlatypusRelease/ Platypus-0.1.5
rm Platypus_0.1.5.tgz; cd ..; echo

# Download transcript database
message "Downloading transcript database"
wget http://www.well.ox.ac.uk/bioinformatics/Software/default_transcripts.zip
unzip default_transcripts.zip
rm default_transcripts.zip; echo

# Build BWA
message "Building BWA"
cd tools/bwa-0.5.10; make; cd ../..; echo 

# Build Stampy
message "Building Stampy"
cd tools/stampy-1.0.14.1
cp makefile makefile.original
sed 's/-Wl//g' makefile.original > makefile
make; cd ../..; echo 

# Build and install Platypus
message "Building and installing Platypus"
cd tools/Platypus-0.1.5; python setup.py build; echo
mkdir $opexdir/Platypus
python setup.py install --prefix $opexdir/Platypus; cd ../..
echo -e "import sys\nsys.path.insert(1, '"$opexdir"/Platypus/lib/python2.7/site-packages')" > usercustomize.py

# Build Pysam 0.7.7
message "Building Pysam"
./install_pysam.sh

mkdir tmp

