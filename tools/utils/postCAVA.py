from __future__ import division
import sys
import os
from subprocess import call

fn=sys.argv[1]
prefn=fn[:-4]+'_pre.vcf'
call(['mv',fn,prefn])
out=open(fn,'w')
for line in open(prefn):
	if line.startswith('#CHROM'):
		cols=line.split('\t')
		x=cols[-1]
		cols[-1]=x[:x.find('_picard.bam')]
		out.write('\t'.join(cols)+'\n')
	else:
		out.write(line)
out.close()
call(['rm',prefn])
	
header=['CHROM','POS','REF','ALT','QUAL','QUALFLAG','FILTER','TR','TC','SAMPLE','GT','TYPE','ENST','GENE','TRINFO','LOC','CSN','CLASS','SO','IMPACT','ALTANN','ALTCLASS','ALTSO']
print '#'+'\t'.join(header)

for line in open(fn):
	line=line.strip()
	if line=='': continue
	if line.startswith('##'): continue
	
	cols=line.split('\t')
	
	if line.startswith('#'): 
		sample=cols[9]
		continue	
	
	chrom=cols[0]
	if chrom.startswith('chr'): chrom=chrom[3:]
	pos=cols[1]
	id=cols[2]
	ref=cols[3]
	alts=cols[4].split(",")
	qual=cols[5]
	filter=cols[6]
	info=cols[7]
	
	infobits=info.split(';')
	infodict=dict()
	for infobit in infobits:
		idx=infobit.find('=')
		if idx!=-1:
			key=infobit[:idx].strip()
			value=infobit[idx+1:].strip()
			infodict[key]=value
			
	ENST_byalt=infodict['ENST'].split(',')
	TYPE_byalt=infodict['TYPE'].split(',')
	GENE_byalt=infodict['GENE'].split(',')
	TRINFO_byalt=infodict['TRINFO'].split(',')
	LOC_byalt=infodict['LOC'].split(',')
	CSN_byalt=infodict['CSN'].split(',')
	CLASS_byalt=infodict['CLASS'].split(',')
	SO_byalt=infodict['SO'].split(',')	
	IMPACT_byalt=infodict['IMPACT'].split(',')
	ALTANN_byalt=infodict['ALTANN'].split(',')
	ALTCLASS_byalt=infodict['ALTCLASS'].split(',')
	ALTSO_byalt=infodict['ALTSO'].split(',')
	
	TRs=infodict['TR'].split(',')
	TC=infodict['TC']
	
	if float(TC)==0: continue
	
	GT=cols[9][:3]
	
	for i in range(len(alts)):
		alt=alts[i]	
		ENST=ENST_byalt[i]
		transcripts=ENST.split(':')
		GENE=GENE_byalt[i]
		GENE_bytrans=GENE.split(':')
		TRINFO=TRINFO_byalt[i]
		TRINFO_bytrans=TRINFO.split(':')
		LOC=LOC_byalt[i]
		LOC_bytrans=LOC.split(':')	
		CSN=CSN_byalt[i]
		CSN_bytrans=CSN.split(':')		
		CLASS=CLASS_byalt[i]
		CLASS_bytrans=CLASS.split(':')
		SO=SO_byalt[i]
		SO_bytrans=SO.split(':')
		IMPACT=IMPACT_byalt[i]
		IMPACT_bytrans=IMPACT.split(':')
		ALTANN=ALTANN_byalt[i]
		ALTANN_bytrans=ALTANN.split(':')
		ALTCLASS=ALTCLASS_byalt[i]
		ALTCLASS_bytrans=ALTCLASS.split(':')
		ALTSO=ALTSO_byalt[i]
		ALTSO_bytrans=ALTSO.split(':')
		
		qualflag=''
		if TYPE_byalt[i]=='Substitution':
			if float(qual)>=100: qualflag='high'
			else: qualflag='low'
		else:
			prop=float(TRs[i])/float(TC)
			if prop>0.2 and filter=='PASS': qualflag='high'
			else: qualflag='low'
			
		for j in range(len(transcripts)):
			
			if CSN_bytrans[j]=='.': continue
			
			out=[chrom,pos,ref,alt,qual,qualflag,filter,TRs[i],TC,sample,GT,TYPE_byalt[i],transcripts[j],GENE_bytrans[j],TRINFO_bytrans[j],LOC_bytrans[j],CSN_bytrans[j],CLASS_bytrans[j],SO_bytrans[j],IMPACT_bytrans[j],ALTANN_bytrans[j],ALTCLASS_bytrans[j],ALTSO_bytrans[j]]
			
			print '\t'.join(out)
		


