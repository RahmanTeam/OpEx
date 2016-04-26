
import json
from collections import OrderedDict

#######################################################################################################################

# Class representing a single Ensembl transcript
class Transcript(object):

	# Constructor
	def __init__(self,line):
		self.exons=[]
		cols=line.split('\t')
		self.ENST=cols[0]
		self.geneSymbol=cols[1]
		self.chrom=cols[2]
		self.strand=int(cols[3])
		self.transcriptStart=int(cols[4])
		self.transcriptEnd=int(cols[5])
		self.codingStart=int(cols[6])
		self.codingStartGenomic=int(cols[7])
		self.codingEndGenomic=int(cols[8])
		# Initializing and adding exons
		for i in range(1,len(cols)-9,2):
			self.exons.append(Exon(int((i+1)/2),int(cols[8+i]),int(cols[9+i])))

	# Getting basic information about the transcript
	def getInfo(self):
		if self.strand==1: ret='+/'
		else: ret='-/'
		cds=0
		for exon in self.exons: cds+=exon.length
		return ret+str(round((self.transcriptEnd-self.transcriptStart+1)/1000,1))+'kb'+'/'+str(len(self.exons))+'/'+str(round(cds/1000,1))+'kb'

	# Getting the full coding sequence of the transcript
	def getCodingSequence(self,reference,variant):
		ret=''
		# Returning reference coding sequence (no variant)
		if variant is None:
			# Iterating through exons and retrieving their reference sequences
			for exon in self.exons:
				if self.strand==1:
					ret+=reference.getReference(self.chrom,exon.start+1,exon.end)
				else:
					ret+=reference.getReference(self.chrom,exon.start+1,exon.end).reverseComplement()
			ret=ret[self.codingStart-1:]
			return ret

		# Returning coding sequence containing a variant
		# Iterating through exons
		for exon in self.exons:
			# CHecking if the variant is located within the exon
			if exon.start<variant.pos<=exon.end:
				if self.strand==1:
					ret+=reference.getReference(self.chrom,exon.start+1,variant.pos-1)+variant.alt+reference.getReference(self.chrom,variant.pos+len(variant.ref),exon.end)
				else:
					temp=Sequence(reference.getReference(self.chrom,exon.start+1,variant.pos-1)+variant.alt+reference.getReference(self.chrom,variant.pos+len(variant.ref),exon.end))
					ret+=temp.reverseComplement()
				continue
			# Adding the reference sequence of the exon, if the variant is not within the exon
			if self.strand==1:
				ret+=reference.getReference(self.chrom,exon.start+1,exon.end)
			else:
				ret+=reference.getReference(self.chrom,exon.start+1,exon.end).reverseComplement()
		# Remove 5' UTR
		ret=ret[self.codingStart-1:]
		return ret

	# Getting the translated protein sequence of the transcript
	def getProteinSequence(self,reference,variant):
		# Translating coding sequence
		ret=Sequence(self.getCodingSequence(reference,variant)).translate(1)
		stop=ret.find('X')
		if stop!=-1: ret=ret[:stop+1]
		return ret

	# Checking if a given position is within the UTR of the transcript
	def isInUTR(self,pos):
		if self.strand==1:
			return (pos<self.codingStartGenomic) or (pos>self.codingEndGenomic)
		else:
			return (pos>self.codingStartGenomic) or (pos<self.codingEndGenomic)

	# Checking if the given variant is outside of the translated region of the transcript
	def isOutsideTranslatedRegion(self,variant):
		if self.strand==1:
			if variant.isInsertion():
				if variant.pos<=self.codingStartGenomic: return True
				if variant.pos-1>=self.codingEndGenomic: return True
				return False
			else:
				if variant.pos+len(variant.ref)-1<self.codingStartGenomic: return True
				if variant.pos>self.codingEndGenomic: return True
				return False
		else:
			if variant.isInsertion():
				if variant.pos<=self.codingEndGenomic: return True
				if variant.pos-1>=self.codingStartGenomic: return True
				return False
			else:
				if variant.pos+len(variant.ref)-1<self.codingEndGenomic: return True
				if variant.pos>self.codingStartGenomic: return True
				return False
  
	# Checking if the given variant is outside of the translated region of the transcript, +/- the first and last 3 bases of the coding sequence
	def isOutsideTranslatedRegionPlus3(self,variant):
		if self.strand==1:
			if variant.isInsertion():
				if variant.pos<=self.codingStartGenomic+3: return True
				if variant.pos-1>=self.codingEndGenomic-3: return True
				return False
			else:
				if variant.pos+len(variant.ref)-1<self.codingStartGenomic+3: return True
				if variant.pos>self.codingEndGenomic-3: return True
				return False
		else:
			if variant.isInsertion():
				if variant.pos<=self.codingEndGenomic+3: return True
				if variant.pos-1>=self.codingStartGenomic-3: return True
				return False
			else:
				if variant.pos+len(variant.ref)-1<self.codingEndGenomic+3: return True
				if variant.pos>self.codingStartGenomic-3: return True
				return False

	# Checking if the given variant overlaps with splicing region
	def isInSplicingRegion(self,variant,ssrange):
		if not self.isOutsideTranslatedRegion(variant):
			for exon in self.exons:
				if variant.overlap(exon.end+1,exon.end+ssrange): return True
				if variant.overlap(exon.start-(ssrange-1),exon.start): return True
			return False
		else:
			return False

	# Checking if the given variant affects an essential splice site
	def isInEssentialSpliceSite(self,variant):
		if not self.isOutsideTranslatedRegion(variant):
			for exon in self.exons:
				if variant.overlap(exon.end+1,exon.end+2): return True
				if variant.overlap(exon.start-1,exon.start): return True
			return False
		else:
			return False

	# Checking if the given variant affects a +5 essential splice site
	def isIn_SS5_Site(self,variant):
		if not self.isOutsideTranslatedRegion(variant):
			if self.strand==1:
				for exon in self.exons:
					if variant.overlap(exon.end+1,exon.end+5):
						if not variant.isSNP():
							if not (variant.pos==exon.end+3 and len(variant.ref)==2 and len(variant.alt)==2): return True
						else:
							if variant.pos==exon.end+5: return True
				return False
			else:
				for exon in self.exons:
					if variant.overlap(exon.start-4,exon.start):
						if not variant.isSNP():
							if not (variant.pos==exon.start-3 and len(variant.ref)==2 and len(variant.alt)==2): return True
						else:
							if variant.pos==exon.start-4: return True
				return False
		else:
			return False

	# Checking if the given variant affects the first or last 3 bases of an exon
	def isInFirstOrLast3BaseOfExon(self,variant):
		if not self.isOutsideTranslatedRegionPlus3(variant):
			for exon in self.exons:
				if variant.overlap(exon.start+1,exon.start+3): return True
				if variant.overlap(exon.end-2,exon.end): return True
			return False
		else:
			return False

	# Checking where a given genomic position is located in the transcript
	def whereIsThisPosition(self,pos):
		if (self.strand==1 and pos<self.codingStartGenomic) or (self.strand==-1 and pos>self.codingStartGenomic): return '5UTR'
		if (self.strand==1 and pos>self.codingEndGenomic) or (self.strand==-1 and pos<self.codingEndGenomic): return '3UTR'
		# Iterating through exons and introns and checking if genomic position is located within
		for exon in self.exons:
			if exon.index>1 and ((self.strand==1 and prevexonend<pos<=exon.start) or (self.strand==-1 and exon.end<pos<=prevexonend)):
				if self.intronLength(exon.index)>5: return 'In'+str(exon.index-1)+'/'+str(exon.index)
				else: return 'fsIn'+str(exon.index-1)+'/'+str(exon.index)
			if exon.start<pos<=exon.end:
				return 'Ex'+str(exon.index)
			prevexonend = exon.end if self.strand==1 else exon.start

	# Checking where a given variant is located in the transcript
	def whereIsThisVariant(self,variant):
		# Getting the locations of both end points of the variant
		if variant.isInsertion():
			first=self.whereIsThisPosition(variant.pos-1)
			second=self.whereIsThisPosition(variant.pos)
		else:
			first=self.whereIsThisPosition(variant.pos)
			second=self.whereIsThisPosition(variant.pos+len(variant.ref)-1)
		if first==second: return first
		if self.strand==1: return first+'-'+second
		else: return second+'-'+first
		

	# Getting the length of an intron, where idx is the index of the succeeding exon
	def intronLength(self,idx):
		for exon in self.exons:
			if exon.index==idx:
				if self.strand==1: return exon.start-prev
				else: return prev-exon.end
			if self.strand==1: prev=exon.end
			else: prev=exon.start

#######################################################################################################################

# Class representing a single exon
class Exon(object):

	# Constructor
	def __init__(self,index,start,end):
		self.index=index
		self.start=start
		self.end=end
		self.length=end-start

#######################################################################################################################       

def getTranscriptCoordinates(enstdb,chrom,pos):
	ret=OrderedDict()
	transcripts=findTranscripts(enstdb,chrom,pos)
	for enstid,transcript in transcripts.iteritems():
		x,y = transformToCSNCoordinate(pos,transcript)
		transcoord = 'c.'+str(x)
		if y!=0:
			if y>0: transcoord+='+'+str(y)
			else: transcoord+=str(y)
		ret[transcript]=transcoord
	return ret

def findTranscripts(enstdb,chrom,pos):  
	ret=OrderedDict()
	
	if chrom in enstdb.contigs:
		goodchrom=chrom
	else:
		if 'chr'+chrom in enstdb.contigs:
			goodchrom='chr'+chrom
		else:
			if chrom.startswith('chr') and chrom[3:] in enstdb.contigs:
				goodchrom=chrom[3:]
			else:
				return ret

	# Checking both end points of the variant
	reg=goodchrom+':'+str(pos)+'-'+str(pos)
	hits=enstdb.fetch(region=reg)

	for line in hits:
		transcript=Transcript(line)
		if not (transcript.transcriptStart+1<=pos<=transcript.transcriptEnd): continue
		ret[transcript.ENST]=transcript
	  
	return ret
		
# Transforming a genomic position ot CSN coordinate
def transformToCSNCoordinate(pos,transcript):
	prevExonEnd=99999999
	# Checking if genomic position is within translated region
	if not transcript.isInUTR(pos):
		sumOfExonLengths = -transcript.codingStart+1
		# Iterating through exons
		for i in range(len(transcript.exons)):
			exon=transcript.exons[i]
			if i>0:
				if transcript.strand==1:
					# Checking if genomic position is within intron
					if prevExonEnd<pos<exon.start+1:
						if pos<=(exon.start+1-prevExonEnd)/2+prevExonEnd:
							x,y = transformToCSNCoordinate(prevExonEnd,transcript)
							return x,pos-prevExonEnd
						else:
							x,y = transformToCSNCoordinate(exon.start+1,transcript)
							return x,pos-exon.start-1
				else:
					# Checking if genomic position is within intron
					if exon.end<pos<prevExonEnd:
						if pos>=(prevExonEnd-exon.end+1)/2+exon.end:
							x,y = transformToCSNCoordinate(prevExonEnd,transcript)
							return x,prevExonEnd-pos
						else:
							x,y = transformToCSNCoordinate(exon.end,transcript)
							return x,exon.end-pos
			# Checking if genomic position is within exon
			if exon.start+1<=pos<=exon.end:
				if transcript.strand==1: return sumOfExonLengths+pos-exon.start,0
				else: return sumOfExonLengths+exon.end-pos+1,0
			# Calculating sum of exon lengths up to this point
			sumOfExonLengths+=exon.length
			if transcript.strand==1: prevExonEnd=exon.end
			else: prevExonEnd=exon.start+1
	# If genomic position is within UTR
	else:
		if transcript.strand==1:
			if pos<transcript.codingStartGenomic: return pos-transcript.codingStartGenomic,0
			if pos>transcript.codingEndGenomic: return '*'+str(pos-transcript.codingEndGenomic),0
		else:
			if pos>transcript.codingStartGenomic: return transcript.codingStartGenomic-pos,0
			if pos<transcript.codingEndGenomic: return '*'+str(transcript.codingEndGenomic-pos),0
