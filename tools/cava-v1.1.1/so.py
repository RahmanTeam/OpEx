#!/usr/bin/python

"""
@author: Marton Munz
Last updated: 18.02.2015
"""

# Sequence Ontology (SO)
#######################################################################################################################

from basics import Sequence
from basics import Variant

#######################################################################################################################

def getSequenceOntologyAnnotation(variant,transcript,protein,mutprotein):
 	where=transcript.whereIsThisVariant(variant)
    
 	if '5UTR' in where: return '5_prime_UTR_variant'
 	if '3UTR' in where: return '3_prime_UTR_variant'
	
 	if '-' in where:
		first=where[:where.find('-')]
		if first.startswith('In'): return 'splice_acceptor_variant'
		if first.startswith('Ex'): return 'splice_donor_variant'
		if first.startswith('fsIn'): return 'splice_acceptor_variant'

	
	if 'In' in where: 
		
		
		if isInSpliceDonor(transcript,variant): return 'splice_donor_variant'
		if isInSpliceAcceptor(transcript,variant): return 'splice_acceptor_variant'
		
		if transcript.intronLength(int(where[where.find('/')+1:]))>=9:
			if transcript.isIn_SS5_Site(variant): return 'splice_donor_5th_base_variant'
	
	out=[]
	
	if isInSplicingRegion(transcript,variant): 
		if where.startswith('In') or where.startswith('fsIn'): return 'intron_variant|splice_region_variant'
		else: out.append('splice_region_variant')
		
	if where.startswith('In') or where.startswith('fsIn'): return 'intron_variant'
    
	if variant.isInFrame():
		if variant.isDeletion(): out.append('inframe_deletion')
		if variant.isInsertion(): out.append('inframe_insertion')
	else:
		out.append('frameshift_variant')
		return '|'.join(out)

	if protein==mutprotein:
		out.append('synonymous_variant')
		return '|'.join(out)
    
	if protein[0]!=mutprotein[0]: 
		out.append('initiator_codon_variant')
		return '|'.join(out)
    
	if (not protein==mutprotein) and len(protein)==len(mutprotein): out.append('missense_variant')  

	while len(protein)>0 and len(mutprotein)>0: 
		if protein[0]==mutprotein[0]:
			protein = protein[1:]
			mutprotein = mutprotein[1:]
		else: break

	if protein[0]=='X' and mutprotein[0]!='X': out.append('stop_lost')  

	while len(protein)>0 and len(mutprotein)>0:
		if protein[-1]==mutprotein[-1]:
			protein = protein[:-1]
			mutprotein = mutprotein[:-1]
		else: break


	if 'X' in mutprotein: out.append('stop_gained')  
	
	
	if ('stop_gained' in out) or ('stop_lost' in out): 
		if 'missense_variant' in out: out.remove('missense_variant')
    
	return '|'.join(out)
        
def isInSpliceDonor(transcript,variant):
   
    if transcript.strand==1:
        for exon in transcript.exons:
            isLastExon=(exon.index==len(transcript.exons))
            if not isLastExon and variant.overlap(exon.end+1,exon.end+2): return True
        return False
    else:
        for exon in transcript.exons:
            isLastExon=(exon.index==len(transcript.exons))
            if not isLastExon and variant.overlap(exon.start-1,exon.start): return True
        return False                
        
def isInSpliceAcceptor(transcript,variant):
    if transcript.strand==1:
        for exon in transcript.exons:
            isFirstExon=(exon.index==1)
            if not isFirstExon and variant.overlap(exon.start-1,exon.start): return True			
        return False
    else:
        for exon in transcript.exons:
            isFirstExon=(exon.index==1) 
            if not isFirstExon and variant.overlap(exon.end+1,exon.end+2): return True
				
        return False

def isInSplicingRegion(transcript,variant):
    if transcript.strand==1:
        for exon in transcript.exons:
			isFirstExon=(exon.index==1)
			isLastExon=(exon.index==len(transcript.exons))    
			if not isLastExon and variant.overlap(exon.end-2,exon.end+8): return True
			if not isFirstExon and variant.overlap(exon.start-7,exon.start+3): return True
        return False
    else:
		for exon in transcript.exons:
			isFirstExon=(exon.index==1)
			isLastExon=(exon.index==len(transcript.exons))
			if not isLastExon and variant.overlap(exon.start-7,exon.start+3): return True
			if not isFirstExon and variant.overlap(exon.end-2,exon.end+8): return True
		return False
	
	
	


