#!/usr/bin/python

"""
@author: Marton Munz
Last updated: 18.02.2015
"""


# Database preparation tool (dbprep)
#######################################################################################################################


from __future__ import division
import os
import sys
import gzip
import datetime
from optparse import OptionParser
from operator import itemgetter

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + '/../../pysamdir')
import pysam

import urllib 

#######################################################################################################################

# Create compressed output file
def createFile(options):
    if not options.ensembl is None:
        n=generateTranscriptDB(options)
        print '\nA total of '+str(n)+' transcripts have been retrieved.\n'
    else:
        generateSNPDB(options)

# Generate transcript database file 
def generateTranscriptDB(options):
    data=dict()
	
    print 'Ensembl version:  '+options.ensembl
    print 'Reference genome: '+options.genome
	
    if not options.input==None:
        total = numOfRecords(options.input)
        print '\n'+str(total)+' transcript IDs read from '+options.input+'.\n'
    else:
        print '\nAll transcripts will be included.\n'		
	
    sys.stdout.write('Downloading Ensembl database... ')
    sys.stdout.flush()
    url='ftp://ftp.ensembl.org/pub/release-'+options.ensembl+'/gtf/homo_sapiens/Homo_sapiens.'+options.genome+'.'+options.ensembl+'.gtf.gz'
    urllib.urlretrieve(url, 'ensembl_data.gz')
    sys.stdout.write('OK\n') 
	
    transIDs=set()
    if not options.input==None: transIDs=readTranscriptIDs(options.input)

    outfile=open('temp.txt','w')

    sys.stdout.write('Extracting transcript data... ')
    sys.stdout.flush()	 	  
    first=True
    prevenst=''
    for line in gzip.open('ensembl_data.gz','r'):
        line=line.strip()
        if line.startswith('#'): continue
        cols=line.split('\t')
        v=cols[8].split(';')
        enst=''
        for x in v:
            x=x.strip()
            if x.startswith('transcript_id'):
                s=x[x.find('\"')+1:]
                enst=s[:s.find('\"')]
                break

        if not options.input==None: 	
            if not enst in transIDs: continue
	
        if first:
            transcript=dict()
            transcript['ENST']=enst
		
            gene=''
            for x in v:
                x=x.strip()
                if x.startswith('gene_name'):
                    s=x[x.find('\"')+1:]
                    gene=s[:s.find('\"')]
                    break	
            transcript['GENE']=gene
			
            transcript['CHROM']=cols[0]
            if cols[6]=='+': transcript['STRAND']='1'
            else: transcript['STRAND']='-1'
            transcript['POS']=0
            transcript['POSEND']=0
            transcript['EXONS']=[]
	
        # Finalizing and printing out transcript object
        if not enst==prevenst and not first:
            finalizeAndOutput(outfile,transcript)
            transcript=dict()
            transcript['ENST']=enst
			
            gene=''
            for x in v:
                x=x.strip()
                if x.startswith('gene_name'):
                    s=x[x.find('\"')+1:]
                    gene=s[:s.find('\"')]
                    break	
            transcript['GENE']=gene

            transcript['CHROM']=cols[0]
            if cols[6]=='+': transcript['STRAND']='1'
            else: transcript['STRAND']='-1'
            transcript['POS']=0
            transcript['POSEND']=0
            transcript['EXONS']=[]
			
        if cols[2]=='exon': 
            idx=0
            for x in v:
                x=x.strip()
                if x.startswith('exon_number'):	
                    s=x[x.find('\"')+1:]
                    idx=int(s[:s.find('\"')])-1
                    break
			
            start=int(cols[3])-1
            end=int(cols[4])
            if idx>=len(transcript['EXONS']):
                for _ in range(len(transcript['EXONS']),idx+1): 
                    transcript['EXONS'].append(None)
            transcript['EXONS'][idx]={'START': start, 'END': end}
	
        if cols[2]=='start_codon' and cols[7]=='0': 
            if transcript['STRAND']=='1': transcript['CODING_START']=int(cols[3])	
            else: transcript['CODING_START']=int(cols[4])	
		
        if cols[2]=='stop_codon': 
            if transcript['STRAND']=='1': transcript['CODING_END']=int(cols[4])
            else: transcript['CODING_END']=int(cols[3])
	
        prevenst=enst
        if first: first=False
	
    finalizeAndOutput(outfile,transcript)			
			
    outfile.close()		
	
    counter=0		
    for line in open('temp.txt'):
        if not line.startswith('ENST'): continue
        counter+=1
        line.rstrip()
        record=line.split('\t')
        record[4]=int(record[4])
        if record[2] in data.keys():
            data[record[2]].append(record)
        else:
            data[record[2]]=[]
            data[record[2]].append(record)
			
    sys.stdout.write('OK\n') 
    sys.stdout.write('Sorting transcripts... ') 
    sys.stdout.flush()	 	
    sortedRecords=sortRecords(data,4)
    writeToFile(sortedRecords,options.output)
	
    sys.stdout.write('OK\n') 
    sys.stdout.write('Removing temporary files... ')
    sys.stdout.flush()	 	         
    os.remove('temp.txt')   
    os.remove('ensembl_data.gz')   
    sys.stdout.write('OK\n') 
	
    return len(sortedRecords)
	
# Generate SNP database file     
def generateSNPDB(options):
    outfile=open(options.output,'w')
    
    print 'dbSNP build: '+str(options.dbSNP_version) 
    if not options.input is None:
        query = readRecords(options.input)
        print str(len(query))+' SNPs are to be retrieved!'
    else:
        print 'All SNPs are to be retrieved!'
    print ''
    
    if options.dbfile.endswith('.gz'):
        dbf = gzip.open(options.dbfile, "r")
    else:
        dbf = open(options.dbfile)
    
    print 'Database: '+options.dbfile
    print 'Searching database...\n'
    counter=0
    counterw=0
    wasHit=False
    for line in dbf:
        line=line.strip()
        
        if not line.startswith('#'):
            counter+=1
            if counter%100000==0: 
                if wasHit and not options.input is None: print '--------------------------------'
                print "Already checked "+str(counter)+" records in "+options.dbfile+'.'
                wasHit=False
                
            row=line.split('\t')
            chrom=row[0]
            pos=row[1]
            ID=row[2]   
            ref=row[3]
            alts=row[4].split(',')            
            
            if not len(ref)==1: continue
            if not all(len(x)==1 for x in alts): continue
            if not options.input is None:
                if not ID in query: continue
                    
            info=row[7]
            info=info[info.index('dbSNPBuildID=')+13:]
            if ';' in info: build=info[:info.index(';')]
            else: build=info
            
            if not int(build)<=int(options.dbSNP_version):
                continue      
        
            if not wasHit and not options.input is None: print '--------------------------------'
            if not options.input is None: print 'SNP '+str(counterw+1)+':  '+ID+' retrieved.'
            outfile.write('\t'.join([ID,chrom,pos])+'\n')    
            counterw+=1
            wasHit=True
            
            if not options.input is None:
                if counterw==len(query):
                    print '\nAll SNPs have been found!'
                    break
                                                       
    dbf.close()
    outfile.close()
    print ''
    print str(counterw)+' SNPs written to file.'
    print ''
    
# Use Tabix to index output file     
def indexFile(options):
    filename=options.output
    if not options.ensembl is None:
        sys.stdout.write('Compressing output file... ') 
        sys.stdout.flush()
        pysam.tabix_compress(filename,filename+'.gz',force=True)
        sys.stdout.write('OK\n') 	
        sys.stdout.write('Indexing output file... ') 
        sys.stdout.flush()
        pysam.tabix_index(filename+'.gz', seq_col=2, start_col=4, end_col=5, meta_char='#',force=True)
        sys.stdout.write('OK\n')
    else:
        print 'Compressing file...'
        pysam.tabix_compress(filename,filename+'.gz',force=True)
        print 'Indexing file...'
        pysam.tabix_index(filename+'.gz', seq_col=1, start_col=2, end_col=2, meta_char='#',force=True)
         	
# Sort records in file     
def sortRecords(records,idx):
   ret=[]
   chroms=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','MT','X','Y']
   for i in range(len(chroms)):
       chrom=chroms[i]
       if chrom in records.keys():
            records[chrom]=sorted(records[chrom], key=itemgetter(idx)) 
   for i in range(len(chroms)):
       chrom=chroms[i]
       if chrom in records.keys():
            for record in records[chrom]:
                ret.append(record)   
   return ret

# Write records to file    
def writeToFile(sortedRecords,filename):
    outfile=open(filename,'w')
    for record in sortedRecords:
        s=str(record[0]).rstrip()
        for i in range(1,len(record)): s+='\t'+str(record[i]).rstrip()
        outfile.write(s+'\n')
    outfile.close()         

# Count number of records in file
def numOfRecords(inputfn):
    counter=0
    for line in open(inputfn):
        if not line.strip()=='':
            counter+=1
    return counter

# Read records from file as a list
def readRecords(inputfn):
    ret=[]
    for line in open(inputfn):
        ret.append(line.strip())
    return ret

# Check if options are correct    
def checkOptions(options):
    if not options.ensembl is None:
        if not options.dbSNP_version is None:
            print 'ERROR: options -e and -s cannot be used together.\n'
            quit()
        if not options.dbfile is None:
            print 'ERROR: options -e and -d cannot be used together.\n'
            quit()		
        if not options.input is None and not os.path.isfile(options.input):
            print 'ERROR: input file not found.\n'
            quit()  
					
    if not options.dbSNP_version is None:
        if options.dbfile is None:
            print 'ERROR: the dbSNP database (-d) must be given\n'
            quit()
        if not options.input is None and not os.path.isfile(options.input):
            print 'ERROR: input file not found.\n'
            quit()
        if not os.path.isfile(options.dbfile):
            print 'ERROR: dbSNP database file not found.\n'
            quit()
			
    if options.output is None:
        print 'ERROR: the output file name (-o) must be given.\n'
        quit()
		
# Read transcript IDs from file  
def readTranscriptIDs(inputfn):
    ret=set()
    for line in open(inputfn):
        ret.add(line.strip())
    return ret
	
# Finalize and output transcript data
def finalizeAndOutput(outfile,transcript):

	if not 'CODING_START' in transcript.keys() or not 'CODING_END' in transcript.keys(): return
	
	if transcript['STRAND']=='1':

		transcript['POS']=transcript['EXONS'][0]['START']
		transcript['POSEND']=transcript['EXONS'][len(transcript['EXONS'])-1]['END']

		codingStartRelative=0
		for exondata in transcript['EXONS']:
			if exondata['START']<=transcript['CODING_START']<=exondata['END']:
				codingStartRelative+=transcript['CODING_START']-exondata['START']
				break
			else:
				codingStartRelative+=exondata['END']-exondata['START']
		transcript['CODING_START_RELATIVE']=codingStartRelative
	else:

		transcript['POS']=transcript['EXONS'][len(transcript['EXONS'])-1]['START']
		transcript['POSEND']=transcript['EXONS'][0]['END']

		codingStartRelative=0
		for exondata in transcript['EXONS']:
			if exondata['START']<=transcript['CODING_START']<=exondata['END']:
				codingStartRelative+=exondata['END']-transcript['CODING_START']+1
				break
			else:
				codingStartRelative+=exondata['END']-exondata['START']
		transcript['CODING_START_RELATIVE']=codingStartRelative

	out=transcript['ENST']+'\t'+transcript['GENE']+'\t'+transcript['CHROM']+'\t'+transcript['STRAND']+'\t'+str(transcript['POS'])
	out+='\t'+str(transcript['POSEND'])+'\t'+str(transcript['CODING_START_RELATIVE'])+'\t'+str(transcript['CODING_START'])
	out+='\t'+str(transcript['CODING_END'])
	for exondata in transcript['EXONS']: out+='\t'+str(exondata['START'])+'\t'+str(exondata['END'])
	outfile.write(out+'\n')


#######################################################################################################################

if __name__ == '__main__': 
   
    # Version number
    ver='v1.1.1'
   
    # Command line argument parsing
    descr='dbprep is a simple tool for generating local Ensembl transcript and dbSNP database files being used by CAVA (via the @ensembl and @dbsnp option flags).'
    epilog='\nExample usage 1 (Generating transcript database): \npython path/to/cava/dbprep.py -i trans.txt -e 65 -g GRCh37 -o out \n\nExample usage 2 (Generating SNP database): \npython path/to/cava/dbprep.py -i SNPs.txt -s 137 -d 00-All.vcf.gz -o out \n\n(More details in the documentation: cava-'+ver+'-doc.pdf)\n\n'
    OptionParser.format_epilog = lambda self, formatter: self.epilog
    parser = OptionParser(usage='python path/to/cava/dbprep.py <options>',version=ver,description=descr,epilog=epilog)
    parser.add_option('-i', "--in", default=None, dest='input', action='store', help="Input filename (list of ENST or dbSNP IDs)")
    parser.add_option('-o', "--out", default=None, dest='output', action='store', help="Output filename prefix")
    parser.add_option('-e', "--ens", default=None, dest='ensembl', action='store', help="Ensembl release version") 
    parser.add_option('-g', "--genome", dest='genome', action='store', default='GRCh37', help="Human genome reference version (default: %default)") 
    parser.add_option('-s', "--snp", default=None, dest='dbSNP_version', action='store', help="dbSNP release version")
    parser.add_option('-d', "--dbsnp", default=None, dest='dbfile', action='store', help="Local dbSNP database file (.vcf or .vcf.gz)")
    (options, args) = parser.parse_args()
    checkOptions(options)
    
    # Printing out version information  
    print "\n-----------------------------------------------------------------"
    print 'CAVA '+ver+' dbprep (database preparation script) is now running.'
    print 'Started: ', datetime.datetime.now(),'\n'
    
    # Creating compressed output file
    createFile(options)
        
    # Indexing output file with Tabix
    indexFile(options)
    
    # Removing uncompressed output file
    os.remove(options.output)
    
    # Printing out summary information
    print ''
    print '---------------------'
    print 'Output files created:'
    print '---------------------'
    if options.ensembl is None: info='SNP database'
    else: info='transcript database'
    print options.output+'.gz ('+info+')'
    print options.output+'.gz.tbi (index file)' 
    print ''
    print 'CAVA dbprep successfully finished: ', datetime.datetime.now()
    print "-----------------------------------------------------------------\n"

#######################################################################################################################
