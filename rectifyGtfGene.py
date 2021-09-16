from collections import defaultdict
trans2gene=defaultdict(str)
with open('../mergeSqantiQc/total_classification.txt','r') as refF:
	line=refF.readline()
	while True:
		line=refF.readline()
		if line=="":
			break
		line=line.split('\n')[0]
		line=line.split('\t')
		trans=line[0]
		gene=line[6]
		#if "PS" in gene and len(gene)==11:
		trans2gene[trans]=gene
with open('../mergeSqantiQc/total.merge_corrected.gtf','r') as tarF, \
     open('../mergeSqantiQc/total.asana.gtf','w') as of:
	while True:
		line=tarF.readline()
		if line == "":
			break
		line=line.split('\n')[0]
		line=line.split('\t')
		if 'scaff' in line[0]:
			scaId=int(line[0].split('_')[1])
			scaId+=12
			line[0]="chr"+str(scaId)
		info=line[8].split(';')
		transInfo=info[0]
		trans=transInfo.split('"')[1]
		if trans2gene[trans] =='':
			print(trans)
			continue
		for i in range(8):
			print(line[i],end='\t',file=of)
		print(transInfo,end= ';',file=of)
		print(" gene_id \"%s\"; gene_name \"%s\";" %(trans2gene[trans],trans2gene[trans]),file=of)
