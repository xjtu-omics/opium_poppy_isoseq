from collections import defaultdict
import pandas as pd
import re
from pyfaidx import Fasta
import sys
sys.path.append('/data/home/xutun/mySrc')
from ttlib.basicInfo import bainfo
from ttlib.stringBase import getGc

dd = '/data/home/xutun/poppy/modifyPaper'
gene2transF = f'{dd}/gene2trans.tab'
gene2gcF = f'{dd}/gene2gc.tab'
gene2fpkmF = f'{dd}/gene2fpkm.tab'
gene2annotStatusF = f'{dd}/gene2AnnotStatus.tab'
gene2exonnF = f'{dd}/gene2exonn.tab'
gene2taxidF = f'{dd}/prot.5000.gene.classF'
geneHaveProtF = f'{dd}/geneHaveProt.tab'
classF = f'{dd}/total_classification.txt'
gene2lenF = f'{dd}/gene2len.tab'
gene2exonnF = f'{dd}/gene2exonn.tab'
gene2exonLenF = f'{dd}/gene2exonLen.tab'
gene2intronLenF = f'{dd}/gene2intronLen.tab'
gene2exon5GcDiffF = f'{dd}/gene2exon5Diff.tab'
gene2exon3GcDiffF = f'{dd}/gene2exon3Diff.tab'
gene2blastpHitsF = f'{dd}/gene2BlastpHitsNum.tab'
gene2pfamDomainNumF = f'{dd}/gene2pfamDomainNum.tab'
gene2GOListF = f'{dd}/gene2GOList.tab'
protGeneListF = f'{dd}/totGene.list'
transProt2lenF = f'{dd}/transProt2len.tab'
toMergeGeneF = f'{dd}/toMergeGene.tab'
toSplitGeneF = f'{dd}/splitGene.list.txt'
totGff = f'{dd}/total.merge_corrected.gff'
corGene2statusF = f'{dd}/corGene2status.tab'
corGene2transF = f'{dd}/corGene2trans.tab'
corGene2lenF = f'{dd}/corGene2len.tab'
corGene2exonnF = f'{dd}/corGene2exonn.tab'
corGene2GOListF = f'{dd}/corGene2GOList.tab'
corGene2pfamDomainNumF = f'{dd}/corGene2pfamDomainNum.tab'
corGene2blastpHitsF = f'{dd}/corGene2blastpHitsNum.tab'
corGeneFpkmF = f'{dd}/cor.tAs.stringtie.fpkm'
annotGtf = f'{dd}/ref/poppy_v6.final_revised.gtf'
isoGtf = f'{dd}/ref/total.merge_corrected.gtf'
refFastaF = f'{dd}/ref/poppy_ref.fa'
corTransInterRepeatF = f'{dd}/corTransInterRepeat.tab'
corGene2gcF = f'{dd}/corGene2gc.tab'
isoTrans2taxidF = f'{dd}/prot.5000.classF'
annotTrans2taxidF = f'{dd}/protAnnot.5000.classF'
corGene2exonLenF = f'{dd}/corGene2exonLen.tab'
corGene2intronLenF = f'{dd}/corGene2intronLen.tab'
corGene2exon5GcDiffF = f'{dd}/corGene2exon5Diff.tab'
corGene2exon3GcDiffF = f'{dd}/corGene2exon3Diff.tab'
corGene2fpkmF = 'data/corGene2fpkm.tab'

def getCorGene2Exon53GcDiffList():
    gene2exon5GcDiffList = defaultdict(list)
    gene2exon3GcDiffList = defaultdict(list)
    for line in open(corGene2exon5GcDiffF):
        line = line.split('\t')
        gene2exon5GcDiffList[line[0]].append(float(line[1]))
    for line in open(corGene2exon3GcDiffF):
        line = line.split('\t')
        gene2exon3GcDiffList[line[0]].append(float(line[1]))
    return gene2exon5GcDiffList,gene2exon3GcDiffList

def getCorGene2intronLenList():
    gene2intronLenList = defaultdict(list)
    for line in open(corGene2intronLenF):
        line = line.strip().split('\t')
        gene2intronLenList[line[0]].append(int(line[1]))
    return gene2intronLenList

def getCorGene2exonLenList():
    gene2exonLenList = defaultdict(list)
    for line in open(corGene2exonLenF):
        line = line.strip().split('\t')
        gene2exonLenList[line[0]].append(int(line[1]))
    return gene2exonLenList

def getCorGene2taxid():
    gene2taxid = defaultdict(int)
    gene2trans = getCorGene2trans()
    trans2taxid = getTotTrans2taxid()
    for gene in gene2trans:
        trans = gene2trans[gene]
        taxid = trans2taxid[trans]
        gene2taxid[gene] = taxid
    return gene2taxid

def getTotTrans2taxid():
    trans2taxid = defaultdict(int)
    for line in open(isoTrans2taxidF):
        line = line.strip().split('\t')
        trans2taxid[line[0]] = int(line[1])
    for line in open(annotTrans2taxidF):
        line = line.strip().split('\t')
        trans2taxid[line[0]] = int(line[1])
    return trans2taxid

def getCorGene2gc():
    gene2gc = defaultdict(float)
    for line in open(corGene2gcF):
        line = line.strip().split('\t')
        gene2gc[line[0]] = float(line[1])
    return gene2gc

def getRepeatInterTransSet():
    transSet = defaultdict(int)
    for line in open(corTransInterRepeatF):
        line = line.strip().split('\t')
        transSet[line[4]] = 1
    return transSet

def generateValidGffFromTransSet(gene2trans,outGtf):
    transGeneSet = defaultdict(str)
    for gene in gene2trans:
        transGeneSet[gene2trans[gene]] = gene

    trans2bainfoList = defaultdict(list)
    ref = Fasta(refFastaF)

    cutP = 0.2

    outList = []
    for line in open(isoGtf):
        cpLine = line.strip()
        line = line.strip().split('\t')
        trans = line[8].split('"')[1]
        if trans not in transGeneSet:
            continue
        chr = line[0]
        st = int(line[3])
        ed = int(line[4])
        trans2bainfoList[trans].append(bainfo(chr, st, ed, _trans=cpLine))

    for trans in trans2bainfoList:
        tmpList = trans2bainfoList[trans]
        tmpList.sort()
        num = len(tmpList)
        stP = -1
        for i in range(num):
            info = tmpList[i]
            seq = ref[info.chr][info.st - 1:info.ed]
            gc = getGc(seq)
            if gc >= cutP:
                stP = i
                break
        if stP == -1:
            print(f'{trans} st wrong!')
            continue
        edP = -1
        i = num - 1
        while i >= 0:
            info = tmpList[i]
            seq = ref[info.chr][info.st - 1:info.ed]
            gc = getGc(seq)
            if gc >= cutP:
                edP = i
                break
            i = i - 1
        if edP == -1:
            print(f'{trans} ed wrong!')
            continue
        if edP < stP:
            print(f'{trans} overlap!')
            continue
        for i in range(stP, edP + 1):
            outList.append(tmpList[i].trans)

    with open(outGtf, 'w') as of:
        for line in outList:
            print(line, file=of)
        for line in open(annotGtf):
            if line[0] == '#':
                continue
            cpLine = line.strip()
            line = line.strip().split('\t')
            if line[2] != 'exon':
                continue
            info = line[8]
            gene = re.findall('gene_name\s"(.*)?"', info)[0]
            if gene not in transGeneSet:
                continue
            print(cpLine, file=of)

def getAllTrans2coding():
    trans2coding = defaultdict(str)
    classD = pd.read_table(classF,sep='\t')
    for ind,row in classD.iterrows():
        trans2coding[row['isoform']] = row['coding']
    return trans2coding

def getCorGene2fpkm():
    gene2fpkm = defaultdict(float)
    for line in open(corGene2fpkmF):
        line = line.strip().split('\t')
        gene2fpkm[line[0]] = float(line[1])
    return gene2fpkm

def getCorGene2blastpHitsNum():
    gene2blastpHitsNum = defaultdict(int)
    for line in open(corGene2blastpHitsF):
        line = line.strip().split('\t')
        gene2blastpHitsNum[line[0]] = int(line[1])
    gene2status = getCorGene2status()
    for gene in gene2status:
        if gene not in gene2blastpHitsNum:
            gene2blastpHitsNum[gene] = 0
    return gene2blastpHitsNum

def getCorGene2pfamDomainNum():
    gene2pfamDomainNum = defaultdict(int)
    for line in open(corGene2pfamDomainNumF):
        line = line.strip().split('\t')
        gene2pfamDomainNum[line[0]] = int(line[1])
    gene2status = getCorGene2status()
    for gene in gene2status:
        if gene not in gene2pfamDomainNum:
            gene2pfamDomainNum[gene] = 0
    return gene2pfamDomainNum


def getCorGene2GOList():
    gene2GOList = defaultdict(list)
    for line in open(corGene2GOListF):
        line = line.strip().split('\t')
        gene2GOList[line[0]].append(line[1])
    gene2status = getCorGene2status()
    for gene in gene2status:
        if gene not in gene2GOList:
            gene2GOList[gene] = []
    return gene2GOList

def getCorGene2exonn():
    gene2exonn = defaultdict(int)
    for line in open(corGene2exonnF):
        line = line.strip().split('\t')
        gene2exonn[line[0]] = int(line[1])
    return gene2exonn

def getCorGene2len():
    gene2len = defaultdict(int)
    for line in open(corGene2lenF):
        line = line.strip().split('\t')
        gene2len[line[0]] = int(line[1])
    return gene2len

def getCorGene2trans():
    gene2trans = defaultdict(str)
    for line in open(corGene2transF):
        line = line.strip().split('\t')
        gene2trans[line[0]] = line[1]
    return gene2trans

def getCorGene2status():
    gene2status = defaultdict(str)
    for line in open(corGene2statusF):
        line = line.strip().split('\t')
        gene2status[line[0]] = line[1]
    return gene2status

def isMerge(gene):
    gene = gene.split('_')
    if len(gene)>=2:
        if len(gene[0])==11 and len(gene[1])==11:
            if 'PS' in gene[0] and 'PS' in gene[1]:
                return True
    return False

def getToMergeGene2trans():
    mergeGene2trans = defaultdict(str)
    for line in open(toMergeGeneF):
        line = line.split('\t')
        mergeGene2trans[line[0]] = line[1]
    return mergeGene2trans

def getToSplitGene2transList():
    splitGene2transList = defaultdict(list)
    for line in open(toSplitGeneF):
        line = line.strip().split('\t')
        gene = line[0]
        for i in range(1,len(line)):
            splitGene2transList[gene].append(line[i])
    return splitGene2transList

def getTransProt2len():
    transProt2len = defaultdict(int)
    for line in open(transProt2lenF):
        line = line.strip().split('\t')
        transProt2len[line[0]] = int(line[1])
    return transProt2len

def getGeneSet():
    geneSet = defaultdict(int)
    for line in open(protGeneListF):
        line = line.strip()
        geneSet[line] = 1
    return geneSet

def getGene2trans():
    gene2trans = defaultdict(str)
    for line in open(gene2transF):
        line = line.strip().split('\t')
        gene2trans[line[0]] = line[1]
    return gene2trans

def getGene2gc():
    gene2gc = defaultdict(float)
    for line in open(gene2gcF):
        line = line.strip().split('\t')
        gene2gc[line[0]] = float(line[1])
    return gene2gc

def getGene2fpkm():
    gene2fpkm = defaultdict(float)
    for line in open(gene2fpkmF):
        line = line.strip().split('\t')
        gene2fpkm[line[0]] = float(line[1])
    return gene2fpkm

def getGene2annotStatus():
    gene2status = defaultdict(str)
    for line in open(gene2annotStatusF):
        line = line.strip().split('\t')
        gene2status[line[0]] = line[1]
    return gene2status

def getGene2exonn():
    gene2exonn = defaultdict(str)
    for line in open(gene2exonnF):
        line = line.strip().split('\t')
        gene2exonn[line[0]] = int(line[1])
    return gene2exonn

def getGene2taxid():
    gene2taxid = defaultdict(str)
    for line in open(gene2taxidF):
        line = line.strip().split('\t')
        gene2taxid[line[0]] = line[1]
    return gene2taxid

def getGeneHaveProt():
    geneHaveProt = defaultdict(int)
    for line in open(geneHaveProtF):
        line = line.strip().split('\t')
        if line[1] == '1':
            geneHaveProt[line[0]] = 1
    return geneHaveProt

def getTotTrans2gene():
    trans2gene = defaultdict(str)
    classD = pd.read_table(classF,sep='\t')
    for ind, row in classD.iterrows():
        trans2gene[row['isoform']] = row['associated_gene']
    return trans2gene

def getGene2len():
    gene2len = defaultdict(str)
    for line in open(gene2lenF):
        line = line.strip().split('\t')
        gene2len[line[0]] = int(line[1])
    return gene2len

def getGene2exonLenList():
    gene2exonLenList = defaultdict(list)
    for line in open(gene2exonLenF):
        line = line.strip().split('\t')
        gene2exonLenList[line[0]].append(int(line[1]))
    return gene2exonLenList

def getGene2intronLenList():
    gene2intronLenList = defaultdict(list)
    for line in open(gene2intronLenF):
        line = line.strip().split('\t')
        gene2intronLenList[line[0]].append(int(line[1]))
    return gene2intronLenList

def getGene2Exon53GcDiffList():
    gene2exon5GcDiffList = defaultdict(list)
    gene2exon3GcDiffList = defaultdict(list)
    for line in open(gene2exon5GcDiffF):
        line = line.split('\t')
        gene2exon5GcDiffList[line[0]].append(float(line[1]))
    for line in open(gene2exon3GcDiffF):
        line = line.split('\t')
        gene2exon3GcDiffList[line[0]].append(float(line[1]))
    return gene2exon5GcDiffList,gene2exon3GcDiffList

def getGene2blastpHitsNum():
    gene2blastpHitsNum = defaultdict(int)
    for line in open(gene2blastpHitsF):
        line = line.strip().split('\t')
        gene2blastpHitsNum[line[0]] = int(line[1])
    gene2status = getGene2annotStatus()
    for gene in gene2status:
        if gene not in gene2blastpHitsNum:
            gene2blastpHitsNum[gene] = 0
    return gene2blastpHitsNum

def getTransGeneSet():
    transGeneSet = defaultdict(str)
    gene2trans = getGene2trans()
    for gene in gene2trans:
        transGeneSet[gene2trans[gene]] = gene
    return transGeneSet

def getCorTransGeneSet():
    transGeneSet = defaultdict(str)
    gene2trans = getCorGene2trans()
    for gene in gene2trans:
        transGeneSet[gene2trans[gene]] = gene
    return transGeneSet

def getGene2pfamDomainNum():
    gene2pfamDomainNum = defaultdict(int)
    for line in open(gene2pfamDomainNumF):
        line = line.strip().split('\t')
        gene2pfamDomainNum[line[0]] = int(line[1])
    gene2status = getGene2annotStatus()
    for gene in gene2status:
        if gene not in gene2pfamDomainNum:
            gene2pfamDomainNum[gene] = 0
    return gene2pfamDomainNum

def getGene2GOList():
    gene2GOList = defaultdict(list)
    for line in open(gene2GOListF):
        line = line.strip().split('\t')
        gene2GOList[line[0]].append(line[1])
    gene2status = getGene2annotStatus()
    for gene in gene2status:
        if gene not in gene2GOList:
            gene2GOList[gene] = []
    return gene2GOList

def getGene2transBainfoList():
    transProt2len = getTransProt2len()
    trans2gene = getTotTrans2gene()
    gene2transBainfoList = defaultdict(list)
    for line in open(totGff):
        if line[0]=='#':
            continue
        ll = line.split('\t')
        if ll[2]!='transcript':
            continue
        trans = re.findall('ID=(.*?);',ll[8])[0]
        if trans not in transProt2len:
            continue
        gene = trans2gene[trans]
        tmpBainfo = bainfo(ll[0],int(ll[3]),int(ll[4]),_trans=trans,_strand=ll[6])
        gene2transBainfoList[gene].append(tmpBainfo)
    return gene2transBainfoList

def getSelectedNuclFasta(gene2trans,outFile):
    annotNuclFasta = f'{dd}/ref/annot.nuc.fasta'
    isoseqNuclFasta = f'{dd}/total.merge_corrected.fasta'
    transGeneSet = defaultdict(str)
    for gene in gene2trans:
        transGeneSet[gene2trans[gene]] = gene
    with open(outFile,'w') as of:
        outFlag = False
        for line in open(annotNuclFasta):
            line = line.strip()
            if line[0] == '>':
                trans = line.split()[0][1:]
                if trans in transGeneSet:
                    outFlag = True
                else:
                    outFlag = False
                if outFlag:
                    gene = transGeneSet[trans]
                    gene = gene.replace('/','__')
                    print(f'>{gene}',file=of)
            elif outFlag:
                print(line,file=of)
        for line in open(isoseqNuclFasta):
            line = line.strip()
            if line[0] == '>':
                trans = line.split()[0][1:]
                if trans in transGeneSet:
                    outFlag = True
                else:
                    outFlag = False
                if outFlag:
                    gene = transGeneSet[trans]
                    gene = gene.replace('/', '__')
                    print(f'>{gene}',file=of)
            elif outFlag:
                print(line, file=of)


























