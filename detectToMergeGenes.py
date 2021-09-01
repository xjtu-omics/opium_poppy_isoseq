import pandas as pd
from collections import defaultdict
from pyfaidx import Faidx
import sys
sys.path.append('/data/home/xutun/mySrc/modifyPoppyPaper')
from getUse import dd,classF

annotDiamondF = f'{dd}/isoseqDiamondAnnot.prot.diamond'
transProtFa = f'{dd}/total.merge_corrected.faa'
annotProtFa = f'{dd}/ref/poppy_v6.proteins.final_revised.fasta'
toMergeGeneF = f'{dd}/toMergeGene.tab'

classD = pd.read_table(classF,sep='\t')

candidateTransSet = defaultdict(list)
trans2annotDiamondSet = defaultdict(lambda:defaultdict(int))
transProtFaHandle = Faidx(transProtFa)
annotProtFaHandle = Faidx(annotProtFa)
resultTrans2protLen = defaultdict(int)


def isMerge(gene):
    gene = gene.split('_')
    if len(gene)>=2:
        if len(gene[0])==11 and len(gene[1])==11:
            if 'PS' in gene[0] and 'PS' in gene[1]:
                return True
    return False

def getTrans2gene():
    trans2gene = defaultdict(int)
    for ind,row in classD.iterrows():
        trans2gene[row['isoform']] = row['associated_gene']
    return trans2gene

def getCandidateCatTransSet():
    transSet = defaultdict(list)
    for ind,row in classD.iterrows():
        if row['coding'] != 'coding':
            continue
        gene = row['associated_gene']
        trans = row['isoform']
        if isMerge(gene):
            transSet[trans] = gene.split('_')
    return transSet

def getTrans2AnnotDiamondSet():
    trans2diamondSet = defaultdict(lambda: defaultdict(int))
    for line in open(annotDiamondF):
        line = line.split('\t')
        trans2diamondSet[line[0]][line[1]] = 1
    return trans2diamondSet

candidateTransSet = getCandidateCatTransSet()
trans2annotDiamondSet = getTrans2AnnotDiamondSet()
trans2gene = getTrans2gene()
relMergeGene2protLen = defaultdict(int)
relMergeGene2trans = defaultdict(str)

for trans in candidateTransSet:
    allInclude = True
    for gene in candidateTransSet[trans]:
        if gene not in trans2annotDiamondSet[trans]:
            allInclude = False
            break
    # if not allInclude:
    #     continue
    allLonger = True
    for gene in candidateTransSet[trans]:
        if transProtFaHandle.index[trans].rlen <= annotProtFaHandle.index[gene].rlen:
            allLonger = False
            break
    # if not allLonger:
    #     continue
    gene = trans2gene[trans]
    if gene not in relMergeGene2protLen or \
       relMergeGene2protLen[gene] < transProtFaHandle.index[trans].rlen:
        relMergeGene2protLen[gene] = transProtFaHandle.index[trans].rlen
        relMergeGene2trans[gene] = trans

with open(toMergeGeneF,'w') as of:
    for gene in relMergeGene2trans:
        print(gene,relMergeGene2trans[gene],sep='\t',file=of)












































