from collections import defaultdict
import pandas as pd
import re
from pyfaidx import Fasta

#输入文件名字,返回文件行数
def countFileLine(fileName):
    return sum(1 for _ in open(fileName))

#提取gtf文件中所有标签为'gene_name'的集合
def gtfGetGeneNameSet(gtfFileName):
    geneNameSet = defaultdict(int)
    for line in open(gtfFileName,'r'):
        line = line.split('\t')
        info = line[8].split(';')
        for item in info:
            if 'gene_name' in item:
                tmpName = item.split('"')[1]
                geneNameSet[tmpName] = 1
    return geneNameSet

#Read list file, seprated by , inside line
def readList(geneListF):
    geneList = []
    for line in open(geneListF,'r'):
        line = line.strip().split(',')
        for ite in line:
            geneList.append(ite)
    return geneList

def writeList(geneList,geneListF):
    with open(geneListF,'w') as of:
        for ite in geneList:
            print(ite,file=of)

def getId2lenFromFai(faiF):
    id2len = defaultdict(int)
    with open(faiF, 'r') as f:
        line = f.readline().strip()
        while line != '':
            line = line.split('\t')
            trans = line[0]
            length = int(line[1])
            id2len[trans] = length
            line = f.readline().strip()
    return id2len

def readMd5OnPd(md5F):
    md5d = pd.read_table(md5F, sep='\s+', header=None, index_col=None)
    md5d.columns = ['md5', 'file']
    return md5d

def readKoTxt(kotxtf):
    gene2keggId = defaultdict(list)
    for line in open(kotxtf):
        line = line.strip().split('\t')
        if len(line) > 1:
            gene2keggId[ line[0] ].append(line[1])
    return gene2keggId

def readWholeCollinearity(colF):
    nowBlockId = 0
    geneA = []
    geneB = []
    blockId = []
    for line in open(colF):
        line = line.strip()
        tmp = re.findall('##\sAlignment\s([0-9])+',line)
        if len(tmp) > 0:
            nowBlockId = int(tmp[0])
        if '#' not in line:
            line = line.split('\t')
            geneA.append(line[1])
            geneB.append(line[2])
            blockId.append(nowBlockId)
    return geneA,geneB,blockId

def fastaGetNameSet(fastaF):
    fastaPandle = Fasta(fastaF)
    nameSet = defaultdict(int)
    for name in fastaPandle.keys():
        nameSet[name] = 1
    return nameSet

def readBuscoFullTsv(tsvF):
    completeBuscoId2GeneList = defaultdict(list)
    duplicatedBuscoId2GeneList = defaultdict(list)
    fragmentedBuscoId2GeneList = defaultdict(list)
    for line in open(tsvF):
        line = line.strip()
        if line[0]=='#':
            continue
        line = line.split('\t')
        buscoId = line[0]
        typ = line[1]
        if typ=='Missing':
            continue
        gene = line[2].split(':')[0]
        if line[1]=='Complete':
            completeBuscoId2GeneList[buscoId].append(gene)
        elif line[1]=='Duplicated':
            duplicatedBuscoId2GeneList[buscoId].append(gene)
        elif line[1]=='Fragmented':
            fragmentedBuscoId2GeneList[buscoId].append(gene)
    return duplicatedBuscoId2GeneList,completeBuscoId2GeneList,duplicatedBuscoId2GeneList