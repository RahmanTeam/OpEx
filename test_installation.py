#!/usr/bin/env python

"""
Last updated: 15.06.2016
"""

from subprocess import call
import os
import sys
import shutil
import filecmp
import os.path


def isReferenceAdded(scriptdir):
    for line in open(scriptdir + '/config.txt'):
        line = line.strip()
        if line == '': continue
        if line.startswith('#'): continue
        if not '=' in line: continue
        [key, value] = line.split('=')
        key = key.strip()
        value = value.strip()
        if key == 'REFERENCE' and value not in ['','.']: return True
    return False

def compareVCF(vcf_x, vcf_y):
    content_x = []
    content_y = []
    for line in open(vcf_x):
        line = line.strip()
        if not line.startswith('#'): content_x.append(line)
    for line in open(vcf_y):
        line = line.strip()
        if not line.startswith('#'): content_y.append(line)
    return (content_x == content_y)



scriptdir = os.path.dirname(os.path.realpath(__file__))
os.chdir(scriptdir)

if not os.path.isfile(scriptdir + '/config.txt'):
    print '\nOpEx not yet installed. Please install it first (See Section 3 in the Documentation).\n'
    quit()

if not isReferenceAdded(scriptdir):
    print '\nReference genome not yet added to the configuration file. You must have performed Quick or Manual installation.'
    print 'Please run Full Installation (see Section 3.2 in the Documentation) or add the reference genome manually (Section 7.4).\n'
    quit()

print '\n'+'='*80
sys.stdout.write('Checking OpEx installation ... ')
sys.stdout.flush()

if os.path.isdir('_testinstall'): shutil.rmtree('_testinstall')

FNULL = open(os.devnull, 'w')

call(['mkdir','_testinstall'])
call(['cp','test/test_R1.fastq.gz','_testinstall'])
call(['cp','test/test_R2.fastq.gz','_testinstall'])
call(['cp','test/test.bed','_testinstall'])

os.chdir('_testinstall')
call(['python','../opex.py','-i','test_R1.fastq.gz,test_R2.fastq.gz','-b','test.bed','-o','test'], stdout=FNULL, stderr=FNULL)
os.chdir('..')

outputs = { '_testinstall/test_calls.vcf': 'test/output/test_calls.vcf',
            '_testinstall/test_annotated_calls.vcf': 'test/output/test_annotated_calls.vcf',
            '_testinstall/test_annotated_calls.txt': 'test/output/test_annotated_calls.txt',
            '_testinstall/test_coverview_summary.txt': 'test/output/test_coverview_summary.txt',
            '_testinstall/test_coverview_regions.txt': 'test/output/test_coverview_regions.txt',
            '_testinstall/test_coverview_profiles.txt': 'test/output/test_coverview_profiles.txt',
            '_testinstall/test_coverview_poor.txt': 'test/output/test_coverview_poor.txt' }

for output, expected in outputs.iteritems():
    if not os.path.isfile(output) or not os.path.isfile(expected):
        print '- Done.'
        print 'OpEx is not installed correctly.'
        print '='*80+'\n'
        quit()

    if output.endswith('.vcf'): same = compareVCF(output, expected)
    else: same = filecmp.cmp(output, expected)
    if not same:
        print '- Done.'
        print 'OpEx is not installed correctly.'
        print '='*80+'\n'
        quit()

print '- Done.'
print 'OpEx is correctly installed.'
shutil.rmtree('_testinstall')









