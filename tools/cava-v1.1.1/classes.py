#!/usr/bin/python

"""
@author: Marton Munz
Last updated: 18.02.2015
"""


# CLASS annotation
#######################################################################################################################

# Getting CLASS annotation of a given variant
def getClassAnnotation(variant,transcript,protein,mutprotein,loc,ssrange):

	# Variants in UTR
	if '5UTR' in loc: return '5PU'
	if '3UTR' in loc: return '3PU'

	# Variants crossing exon-intron boundaries
	if '-' in loc: return 'ESS'

	# Intronic, splice site and essential splice site variants
	if 'In' in loc:
		if transcript.isInEssentialSpliceSite(variant): return 'ESS'
		if transcript.intronLength(int(loc[loc.find('/')+1:]))>=9:
			if transcript.isIn_SS5_Site(variant): return 'SS5'
		if transcript.isInSplicingRegion(variant,ssrange): return 'SS'
		return 'INT'

	# Checking if variant affect the first or last three bases of an exon
	potSS=transcript.isInFirstOrLast3BaseOfExon(variant)

	protL=len(protein)
	mutprotL=len(mutprotein)

	# Frame-shift coding variants
	if not variant.isInFrame(): return 'FS'

	# Synonymous coding variants
	if protein==mutprotein:
		if potSS: return 'EE'
		return 'SY'
		
	# Variants affecting the initiation amino acid
	if protein[0]!=mutprotein[0]: return 'IM'

	# Stop gain and stop lost variants
	while len(protein)>0 and len(mutprotein)>0:
		if protein[0]==mutprotein[0]:
			protein = protein[1:]
			mutprotein = mutprotein[1:]
		else: break



	if protein[0]=='X' and mutprotein[0]!='X': return 'SL'


	while len(protein)>0 and len(mutprotein)>0:
		if protein[-1]==mutprotein[-1]:
			protein = protein[:-1]
			mutprotein = mutprotein[:-1]
		else: break


	if 'X' in mutprotein: return 'SG'
 	

	# Non-synonymous coding variants
	if protL==mutprotL:
		if potSS: return 'EE'
		return 'NSY'

	# In-frame variants
	if potSS: return 'EE'
	return 'IF'

#######################################################################################################################
