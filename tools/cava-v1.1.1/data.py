#!/usr/bin/python

"""
@author: Marton Munz
Last updated: 18.02.2015
"""


# Classes providing interfaces with annotation databases and the reference genome
#######################################################################################################################


import logging
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + '/../../pysamdir')
import pysam

import basics
import classes
import csn
import so


#######################################################################################################################

# Class representing the Ensembl (transcript) dataset
class Ensembl(object):
    
    # Constructor
    def __init__(self,options):
        self.options=options
        # Openning tabix file representing the Ensembl database
        self.tabixfile=pysam.Tabixfile(options.args['ensembl'])
        self.proteinSeqs=dict()
        self.exonSeqs=dict()

    # Finding trancsripts that overlap with a variant 
    def findTranscripts(self,variant_plus,variant_minus,reference):  
        ret=dict()
        retOUT=dict()
        
        # Checking chromosome name
        goodchrom=variant_plus.chrom
        if not goodchrom in self.tabixfile.contigs:
            goodchrom='chr'+goodchrom
            if not goodchrom in self.tabixfile.contigs: return ret,retOUT
        
        if variant_plus.pos==variant_minus.pos: difference=False
        else: difference=True
        
        # Looking at the variant when aligned on the forward strand
        
        # Defining variant end points      
        if not variant_plus.isInsertion():
            start=variant_plus.pos
            end=variant_plus.pos+len(variant_plus.ref)-1
        else:
            start=variant_plus.pos-1
            end=variant_plus.pos
        
        # Checking both end points of the variant
        reg1=goodchrom+':'+str(start)+'-'+str(start)
        reg2=goodchrom+':'+str(end)+'-'+str(end)
        
        
        if not variant_plus.isSNP():
            hits1=self.tabixfile.fetch(region=reg1)
            hits2=self.tabixfile.fetch(region=reg2) 
            hitdict1=dict()
            hitdict2=dict()
            for line in hits1:
                transcript=basics.Transcript(line,reference)
                if not (transcript.transcriptStart+1<=start<=transcript.transcriptEnd): continue
                if difference and transcript.strand==-1: continue
                hitdict1[transcript.ENST]=transcript
            for line in hits2:
                transcript=basics.Transcript(line,reference)
                if not (transcript.transcriptStart+1<=end<=transcript.transcriptEnd): continue
                if difference and transcript.strand==-1: continue
                hitdict2[transcript.ENST]=transcript
		
            # Find transcripts with which the variant fully or partially overlaps
            for key,transcript in hitdict1.iteritems():
                if key in hitdict2.keys():
                    ret[key]=transcript
                else:
                    if not variant_plus.isInsertion(): retOUT[key]=transcript
        
            if not variant_plus.isInsertion():
                for key,transcript in hitdict2.iteritems():
                    if not key in hitdict1.keys():
                        retOUT[key]=transcript
						
        else:
            hits1=self.tabixfile.fetch(region=reg2)
            for line in hits1:
                transcript=basics.Transcript(line,reference)
                if not (transcript.transcriptStart+1<=end<=transcript.transcriptEnd): continue
                ret[transcript.ENST]=transcript
					        
        # Return transcripts if variant alignment makes no difference
        if not difference: return ret,retOUT
        
        # Looking at the variant when aligned on the reverse strand
        
        # Defining variant end points      
        if not variant_minus.isInsertion():
            start=variant_minus.pos
            end=variant_minus.pos+len(variant_minus.ref)-1
        else:
            start=variant_minus.pos-1
            end=variant_minus.pos
        
        # Checking both end points of the variant
        reg1=goodchrom+':'+str(start)+'-'+str(start)
        reg2=goodchrom+':'+str(end)+'-'+str(end)
        hits1=self.tabixfile.fetch(region=reg1)
        hits2=self.tabixfile.fetch(region=reg2) 
        hitdict1=dict()
        hitdict2=dict()
        for line in hits1:
            transcript=basics.Transcript(line,reference)
            if not (transcript.transcriptStart+1<=start<=transcript.transcriptEnd): continue
            if transcript.strand==1: continue
            hitdict1[transcript.ENST]=transcript
        for line in hits2:
            transcript=basics.Transcript(line,reference)
            if not (transcript.transcriptStart+1<=end<=transcript.transcriptEnd): continue
            if transcript.strand==1: continue
            hitdict2[transcript.ENST]=transcript
          
        # Find transcripts with which the variant fully or partially overlaps
        for key,transcript in hitdict1.iteritems():
            if key in hitdict2.keys():
                ret[key]=transcript
            else:
                 if not variant_minus.isInsertion(): retOUT[key]=transcript
        if not variant_minus.isInsertion():
            for key,transcript in hitdict2.iteritems():
                if not key in hitdict1.keys():
                    retOUT[key]=transcript
                              
        return ret,retOUT 
    
    # Returning True if variant is a duplication overlapping minus SS boundary (and for some other INTR calls)
    def isDupOverlappingMinusSSBoundary(self,csn):
        ssrange=int(self.options.args['ssrange'])
        csn_str=csn.getAsString()
        idx=csn_str.find('_p')
        if not idx==-1: cpart=csn_str[:idx]
        else: cpart=csn_str
        idx=cpart.find('dup')
        if idx==-1: return False
        cpart=cpart[2:idx]
        coords=cpart.split('_')
        firstcoord=coords[0]
        idx=firstcoord.find('-')
        if idx==-1:
            if len(coords)==1: return False 
            secondcoord=coords[1]
            idx=secondcoord.find('-')
            if idx<1: return False
            y=int(secondcoord[idx+1:])
            if y<=ssrange: return True
            else: return False
        if idx==0: return False
        x=int(firstcoord[idx+1:])
        if x>=ssrange: return True
        else: return False
		
    # Returning True if variant is a duplication overlapping minus 8
    def isDupOverlappingMinus8Boundary(self,csn):
        ssrange=8
        csn_str=csn.getAsString()
        idx=csn_str.find('_p')
        if not idx==-1: cpart=csn_str[:idx]
        else: cpart=csn_str
        idx=cpart.find('dup')
        if idx==-1: return False
        cpart=cpart[2:idx]
        coords=cpart.split('_')
        firstcoord=coords[0]
        idx=firstcoord.find('-')
        if idx==-1:
            if len(coords)==1: return False 
            secondcoord=coords[1]
            idx=secondcoord.find('-')
            if idx<1: return False
            y=int(secondcoord[idx+1:])
            if y<=ssrange: return True
            else: return False
        if idx==0: return False
        x=int(firstcoord[idx+1:])
        if x>=ssrange: return True
        else: return False
    
    # Annotating a variant based on Ensembl data  
    def annotate(self,variant,reference,impactdir):
        
		
        # Create left-aligned and right-aligned versions of the variant
        if not variant.isSNP():
            variant_plus=variant.alignOnPlusStrand(reference)
            variant_minus=variant.alignOnMinusStrand(reference)
        else:
            variant_plus=variant
            variant_minus=variant
        
        # Checking if variant alignment makes any difference 
        if variant_plus.pos==variant_minus.pos: difference=False
        else: difference=True
        
        # Initializing annotation strings
        ENSTstring=''
        GENEstring=''
        LOCstring=''
        CSNstring=''
        CLASSstring=''
        ALTFLAGstring=''
        TRINFOstring=''
        ALTANNstring=''
        ALTCLASSstring=''
        SOstring=''
        ALTSOstring=''
        IMPACTstring=''
     
        # Collecting transcripts that overlap with the variant
        transcripts,transcriptsOUT=self.findTranscripts(variant_plus,variant_minus,reference)
        
        # Annotating with transcripts that only partially overlap with the variant
        for ENST,transcript in transcriptsOUT.iteritems():        
            if len(ENSTstring)>0:
                ENSTstring+=':'
                GENEstring+=':'
                TRINFOstring+=':'
                LOCstring+=':'
                CSNstring+=':'
                CLASSstring+=':'
                ALTFLAGstring+=':'
                ALTANNstring+=':'
                ALTCLASSstring+=':'
                SOstring+=':'
                ALTSOstring+=':'
                IMPACTstring+=':'
            ENSTstring+=ENST      
            GENEstring+=transcript.geneSymbol
            TRINFOstring+=transcript.getInfo()
            LOCstring+='OUT'
            CSNstring+='.'
            CLASSstring+='.'
            ALTFLAGstring+='.'
            ALTANNstring+='.'
            ALTCLASSstring+='.'
            SOstring+='.'
            ALTSOstring+='.'
            IMPACTstring+='.'
        
        # Iterating through the list of transcripts       
        for ENST,transcript in transcripts.iteritems():
            
            # Separating annotations by different transcripts with colon     
            if len(ENSTstring)>0:
                ENSTstring+=':'
                GENEstring+=':'
                TRINFOstring+=':'
                LOCstring+=':'
                CSNstring+=':'
                CLASSstring+=':'
                ALTFLAGstring+=':'
                ALTANNstring+=':'
                ALTCLASSstring+=':'
                SOstring+=':'
                ALTSOstring+=':'
                IMPACTstring+=':'
               
			   
 			    
            # Creating the ENST, GENE and TRINFO annotations   
            ENSTstring+=ENST      
            GENEstring+=transcript.geneSymbol
            TRINFOstring+=transcript.getInfo()
            
            # Creating the LOC annotation
            loc_plus=transcript.whereIsThisVariant(variant_plus)
            if difference: loc_minus=transcript.whereIsThisVariant(variant_minus)
            else: loc_minus=loc_plus

			# Creating reference and mutated protein sequence
            notexonic_plus=(('5UTR' in loc_plus) or ('3UTR' in loc_plus) or ('-' in loc_plus) or ('In' in loc_plus))
            if difference: notexonic_minus=(('5UTR' in loc_minus) or ('3UTR' in loc_minus) or ('-' in loc_minus) or ('In' in loc_minus))
            else: notexonic_minus=notexonic_plus
            if notexonic_plus and notexonic_minus: protein=''
            else:
                if not transcript.ENST in self.proteinSeqs.keys():				 
                     protein,exonseqs=transcript.getProteinSequence(reference,None,None)
					 
                     if len(self.proteinSeqs)>5: 
                          self.proteinSeqs=dict()
                          self.exonSeqs=dict()
					 
                     self.proteinSeqs[transcript.ENST]=protein
                     self.exonSeqs[transcript.ENST]=exonseqs
                else:
                     protein=self.proteinSeqs[transcript.ENST]
                     exonseqs=self.exonSeqs[transcript.ENST]											 
			
            if notexonic_plus: mutprotein_plus=''
            else: 
                mutprotein_plus,exonseqs=transcript.getProteinSequence(reference,variant_plus,exonseqs)
            
            if difference:
                if notexonic_minus: mutprotein_minus=''
                else: 
                    mutprotein_minus,exonseqs=transcript.getProteinSequence(reference,variant_minus,exonseqs)
            else:
                mutprotein_minus=mutprotein_plus								
   
            # Creating the CSN annotations both for left and right aligned variant    
            csn_plus=csn.getAnnotation(variant_plus,transcript,reference,protein,mutprotein_plus)
            if difference: csn_minus=csn.getAnnotation(variant_minus,transcript,reference,protein,mutprotein_minus)
            else: csn_minus=csn_plus
            
            so_plus=''
            so_minus=''
            class_plus=''
            class_minus=''
            impact_plus=''
            impact_minus=''
				
            if not impactdir==None or self.options.args['ontology'].upper() in ['CLASS','BOTH']:
	            # Creating the CLASS annotations both for left and right aligned variant    
	            class_plus=classes.getClassAnnotation(variant_plus,transcript,protein,mutprotein_plus,loc_plus,int(self.options.args['ssrange'])) 
	            if difference: class_minus=classes.getClassAnnotation(variant_minus,transcript,protein,mutprotein_minus,loc_minus,int(self.options.args['ssrange']))
	            else: class_minus=class_plus
	            # CLASS is set to INTR for duplications overlapping minus SS boundary
	            if (not loc_plus=='Ex1') and (not '-' in loc_plus) and (not (loc_plus=='5UTR' or loc_plus=='3UTR')):
	                if self.isDupOverlappingMinusSSBoundary(csn_plus): 
	                         class_plus='INT'
	                         if loc_plus.startswith('Ex'): 
	                             N=int(loc_plus[2:])
	                             loc_plus='In'+str(N-1)+'/'+str(N)							 
	            if (not loc_minus=='Ex1') and (not '-' in loc_minus) and (not (loc_minus=='5UTR' or loc_minus=='3UTR')):					 
	                if self.isDupOverlappingMinusSSBoundary(csn_minus): 
	                         class_minus='INT'
	                         if loc_minus.startswith('Ex'): 
	                             N=int(loc_minus[2:])
	                             loc_minus='In'+str(N-1)+'/'+str(N) 
			
			# Determining the IMPACT flag
            if not impactdir==None:
	            if class_plus in impactdir.keys(): impact_plus=impactdir[class_plus]
	            else: impact_plus='None'
	            if class_minus in impactdir.keys(): impact_minus=impactdir[class_minus]
	            else: impact_minus='None'
	        
            if self.options.args['ontology'].upper() in ['SO','BOTH']:
                # Creating the SO annotations both for left and right aligned variant   
                so_plus=so.getSequenceOntologyAnnotation(variant_plus,transcript,protein,mutprotein_plus)            
                if difference: so_minus=so.getSequenceOntologyAnnotation(variant_minus,transcript,protein,mutprotein_minus)
                else: so_minus=so_plus  
                # SO is set to intron_variant for duplications overlapping minus SS boundary				 
                if (not loc_plus=='Ex1') and (not '-' in loc_plus) and (not (loc_plus=='5UTR' or loc_plus=='3UTR')):
                    if self.isDupOverlappingMinus8Boundary(csn_plus): 
                             so_plus='intron_variant'
                             if loc_plus.startswith('Ex'): 
                                 N=int(loc_plus[2:])
                                 loc_plus='In'+str(N-1)+'/'+str(N)							 
                if (not loc_minus=='Ex1') and (not '-' in loc_minus) and (not (loc_minus=='5UTR' or loc_minus=='3UTR')):					 
                    if self.isDupOverlappingMinus8Boundary(csn_minus): 
                             so_minus='intron_variant'
                             if loc_minus.startswith('Ex'): 
                                 N=int(loc_minus[2:])
                                 loc_minus='In'+str(N-1)+'/'+str(N)
								 				 
								 			
            # Deciding which is the correct CSN and CLASS annotation
            if transcript.strand==1: 
                CSNstring+=csn_plus.getAsString()        
                CLASSstring+=class_plus
                ALTANN=csn_minus.getAsString()
                altCLASS=class_minus
                SOstring+=so_plus
                altSO=so_minus
                LOCstring+=loc_plus
                IMPACTstring+=impact_plus
            else:
                CSNstring+=csn_minus.getAsString()
                CLASSstring+=class_minus
                ALTANN=csn_plus.getAsString()
                altCLASS=class_plus
                SOstring+=so_minus
                altSO=so_plus
                LOCstring+=loc_minus
                IMPACTstring+=impact_minus
				
            
            if self.options.args['givealt']:
                
                # Creating the ALTANN annotation
                if not csn_plus.getAsString()==csn_minus.getAsString(): ALTANNstring+=ALTANN
                else: ALTANNstring+='.'
                
                # Creating the ALTCLASS annotation
                if not class_plus==class_minus: ALTCLASSstring+=altCLASS
                else: ALTCLASSstring+='.'
				
                # Creating the ALTSO annotations
                if not so_plus==so_minus: ALTSOstring+=altSO
                else: ALTSOstring+='.' 
			                   
            else:				
                # Creating the ALTFLAG annotation     
		
                if self.options.args['ontology'].upper() == 'CLASS':
                    if not class_plus==class_minus: ALTFLAGstring+='AnnAndClass'
                    else:
                        if not csn_plus.getAsString()==csn_minus.getAsString(): ALTFLAGstring+='AnnNotClass'
                        else: ALTFLAGstring+='None'
					
                if self.options.args['ontology'].upper() == 'SO':
                    if not so_plus==so_minus: ALTFLAGstring+='AnnAndSO'
                    else:
                        if not csn_plus.getAsString()==csn_minus.getAsString(): ALTFLAGstring+='AnnNotSO'
                        else: ALTFLAGstring+='None'
					
                if self.options.args['ontology'].upper() == 'BOTH':
					if not class_plus==class_minus:
						if not so_plus==so_minus: ALTFLAGstring+='AnnAndClassAndSO'
						else: ALTFLAGstring+='AnnAndClassNotSO'
					else:
						if csn_plus.getAsString()==csn_minus.getAsString(): ALTFLAGstring+='None'
						else:
							if not so_plus==so_minus: ALTFLAGstring+='AnnAndSONotClass'
							else: ALTFLAGstring+='AnnNotClassNotSO'
							

        # Adding annotations to the variant      
        variant.addFlag('ENST',ENSTstring)
        variant.addFlag('GENE',GENEstring)
        variant.addFlag('TRINFO',TRINFOstring)
        variant.addFlag('LOC',LOCstring)
        variant.addFlag('CSN',CSNstring)
		
        if self.options.args['ontology'].upper() in ['CLASS','BOTH']: variant.addFlag('CLASS',CLASSstring)
        if self.options.args['ontology'].upper() in ['SO','BOTH']: variant.addFlag('SO',SOstring)
		
        if not impactdir==None: variant.addFlag('IMPACT',IMPACTstring)	
      
        if self.options.args['givealt']:
            variant.addFlag('ALTANN',ALTANNstring)
            if self.options.args['ontology'].upper() in ['CLASS','BOTH']: variant.addFlag('ALTCLASS',ALTCLASSstring)
            if self.options.args['ontology'].upper() in ['SO','BOTH']: variant.addFlag('ALTSO',ALTSOstring)
        else: 
            variant.addFlag('ALTFLAG',ALTFLAGstring)        
        
        return variant
                     
#######################################################################################################################            

# Class representing the dbSNP dataset
class dbSNP(object):
    
    # Constructor
    def __init__(self,options):  
        # Openning tabix file representing the dbSNP database
        self.tabixfile=pysam.Tabixfile(options.args['dbsnp'])

    # Annotating a variant based on dbSNP data  
    def annotate(self,variant):
        # Checking if variant is a SNP at all
        if variant.isSNP():
            # Fetching data from dbSNP database
            goodchrom=variant.chrom
            if not goodchrom in self.tabixfile.contigs:
                goodchrom='chr'+goodchrom
                if not goodchrom in self.tabixfile.contigs:
                    variant.addFlag('DBSNP','')
                    return variant    
            reg=goodchrom+':'+str(variant.pos)+'-'+str(variant.pos)
            lines=self.tabixfile.fetch(region=reg)
            for line in lines:
                cols=line.split('\t')
                # Adding DBSNP annotation to the variant 
                variant.addFlag('DBSNP',cols[0])
            if not 'DBSNP' in variant.flags: variant.addFlag('DBSNP','')  
        else:
            variant.addFlag('DBSNP','')
        return variant

#######################################################################################################################

# Class representing the reference genome dataset
class Reference(object):
    
    # Constructor
    def __init__(self,options):
        # Openning tabix file representing the reference genome
        self.fastafile=pysam.Fastafile(options.args['reference'])
        
    # Retrieving the sequence of a genomic region    
    def getReference(self,chrom,start,end): 
        # Checking if chromosome name exists
        goodchrom=chrom
        if not goodchrom in self.fastafile.references:
            goodchrom='chr'+chrom
            if not goodchrom in self.fastafile.references: 
                if chrom=='MT': 
                    goodchrom='chrM'
                    if not goodchrom in self.fastafile.references: return None   
                else: return None  
                     
        # Fetching data from reference genome
        if end<start: return basics.Sequence('')
        if start<0: start=1
        
        if pysam.__version__ in ['0.7.7','0.7.8','0.8.0']:
            last=self.fastafile.getReferenceLength(goodchrom)
        else:
            last=self.fastafile.get_reference_length(goodchrom)
		
        if end>last: end=last
        seq = self.fastafile.fetch(goodchrom,start-1,end)
        return basics.Sequence(seq.upper())
        
#######################################################################################################################