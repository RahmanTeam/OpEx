#!/bin/bash

message () {
echo; printf '%.0s=' {1..80}; echo
echo "OpEx: $1"
printf '%.0s=' {1..80}; echo; echo
}

# Index genome by BWA
message "Indexing genome by BWA"
tools/bwa-0.5.10/bwa index -a bwtsw $1; echo

# Index genome by Stampy
message "Indexing genome by Stampy"
mkdir index
tools/stampy-1.0.14.1/stampy.py --species=human --assembly=hg19_ncbi37 -G index/ref $1; echo
tools/stampy-1.0.14.1/stampy.py -g index/ref -H index/ref; echo


