#!/usr/bin/env python

"""
Last updated: 15/06/2016
"""

import os
from optparse import OptionParser
from subprocess import call
from collections import OrderedDict


def readConfigFile(scriptdir, configpath):
    ret = OrderedDict()
    if configpath == None:
        fn = scriptdir + '/config.txt'
    else:
        fn = configpath
    for line in open(fn):
        line = line.strip()
        if line == '': continue
        if line.startswith('#'): continue
        if not '=' in line: continue
        [key, value] = line.split('=')
        key = key.strip()
        value = value.strip()
        ret[key] = value
    return ret


def writeConfigFile(scriptdir, params):
    out = open(scriptdir + '/config.txt', 'w')
    for key, value in params.iteritems():
        out.write(key + ' = ' + value + '\n')
    out.close()


def generateFile(params, fnin, fnout):
    with open(fnout, "wt") as fout:
        with open(fnin, "rt") as fin:
            for line in fin:
                for key, value in params.iteritems():
                    line = line.replace('@' + key, value)
                fout.write(line)

def checkInputs(options):
    if options.fastq is None:
        print '\nInput files not specified.\n'
        quit()
    x = options.fastq.split(',')
    if not len(x) == 2:
        print '\nIncorrect format for option --input.\n'
        quit()
    if not x[0].endswith('.fastq.gz') or not x[1].endswith('.fastq.gz'):
        print '\nInput files must have .fastq.gz format.\n'
        quit()
    if options.name is None:
        print '\nOutput file name not specified.\n'
        quit()



##############################################################################################################

scriptdir = os.path.dirname(os.path.realpath(__file__))
workingdir = os.getcwd()

# Version
ver = '1.0.0'

# Command line argument parsing
descr = 'OpEx (Optimised Exome) pipeline ' + ver + '.'
parser = OptionParser(usage='python opex.py <options>', version=ver, description=descr)
parser.add_option('-i', "--input", default=None, dest='fastq', action='store', help="fastq.gz files")
parser.add_option('-o', "--output", default=None, dest='name', action='store', help="Sample name (output prefix)")
parser.add_option('-b', "--bed", default=None, dest='bed', action='store', help="Bed file")
parser.add_option('-r', "--reference", default=None, dest='reference', action='store', help="Reference genome file")
parser.add_option('-t', "--threads", default=1, dest='threads', action='store', help="Number of processes to use")
parser.add_option('-f', "--full", default=False, dest='full', action='store_true',
                  help="Output full CoverView output [default value: %default]")
parser.add_option('-c', "--config", default=None, dest='config', action='store', help="Configuration file")
parser.add_option('-k', "--keep", default=False, dest='keep', action='store_true', help="Keep temporary files")
(options, args) = parser.parse_args()
checkInputs(options)

# Welcome message
print '\n' + '-' * 80
print 'OpEx pipeline version ' + ver
print '-' * 80 + '\n'

# Read configuration file
params = readConfigFile(scriptdir, options.config)

if not 'REFERENCE' in params.keys():
    if options.reference is None:
        print 'Error: no reference genom provided.'
        quit()

    refdir = os.path.abspath(options.reference)

    os.chdir(scriptdir)

    params['REFERENCE'] = refdir
    params['GENOME_INDEX'] = scriptdir + '/index/ref'
    params['HASH'] = scriptdir + '/index/ref'
    writeConfigFile(scriptdir, params)
    print 'Adding default reference genome ...'
    call(['chmod', '+x', './index_genome.sh'])
    call(['./index_genome.sh', refdir])

    cavaconfig = open(scriptdir + '/cava_config.txt', "a")
    cavaconfig.write('\n# Name of reference genome file\n')
    cavaconfig.write('# Possible values: string | Optional: no\n')
    cavaconfig.write('@reference = ' + refdir + '\n')
    cavaconfig.close()

    os.chdir(workingdir)

# Additional params
params['NAME'] = options.name
params['FASTQ1'], params['FASTQ2'] = options.fastq.split(',')
params['OPEXDIR'] = scriptdir

params['MORECV'] = ''
if not options.full:
    params['MORECV'] = params['MORECV'] + '-c ' + scriptdir + '/CoverView_default.json'
else:
    params['MORECV'] = params['MORECV'] + '-c ' + scriptdir + '/CoverView_full.json'
if not options.bed is None: params['MORECV'] = params['MORECV'] + ' -b ' + options.bed
if int(options.threads) > 1: params['MORECV'] = params['MORECV'] + ' -t ' + str(options.threads)

params['MORECAVA'] = ''
if int(options.threads) > 1: params['MORECAVA'] = params['MORECAVA'] + '-t ' + str(options.threads)

if options.keep:
    params['KEEPREMOVE'] = ''
else:
    params['KEEPREMOVE'] = 'rm -r ' + params['NAME'] + '_tmp'

# Genearate Bash script file
scriptfn = params['NAME'] + '_opex_pipeline.sh'
generateFile(params, scriptdir + '/templates/opex_pipeline_template', scriptfn)
call(['chmod', '+x', scriptfn])
print 'Bash script has been successfully generated: ' + scriptfn

# Run Bash script 
print '\nRunning the ' + scriptfn + ' script ... '
call(['./' + scriptfn])

# Goodbye message
print '\n' + '-' * 80
print 'OpEx pipeline finished.'
print '-' * 80 + '\n'
