import sys
import os
from collections import defaultdict
sys.path.append('/data/home/xutun/mySrc')
sys.path.append('/data/home/xutun/mySrc/modifyPoppyPaper')
from ttlib.basicInfo import overlapA, overlapMax, overlapMin, bainfo
from ttlib.ttDataStructure import ttUFS
from getUse import getGene2transBainfoList
print('a')
gene2transBainfo = getGene2transBainfoList()
print('b')
trans2Bainfo = defaultdict(bainfo)
overMaxDiffP = 0.2
overMinSameP = 0.8
overMaxSameP = 0.8
maxTransSpan = 30000
splitGeneSet = defaultdict(int)
spg2trans = defaultdict(list)
#临时添加供临时使用,只能小批量查询
def getTranscriptFlncNum(trans):
    if 'transcript' in trans:
        return 2
        # return int(os.popen('grep "^%s," ../data/merged.polished.cluster_report.csv|wc -l'%(trans)).readline())
    else:
        return 1

for gene in gene2transBainfo:
    print(gene)
    tmpList = gene2transBainfo[gene]
    length = len(tmpList)
    totalAbNum = 0
    chrom = tmpList[0].chr
    breakPosList = list()
    splitF = False
    for i in range(length):
        diffNum = 0
        maxSameNum = 0
        if tmpList[i].ed - tmpList[i].st + 1 >= maxTransSpan:
            continue
        for j in range(i + 1, length):
            if tmpList[j].ed - tmpList[j].st + 1 >= maxTransSpan:
                continue
            if overlapMax(tmpList[i], tmpList[j]) <= overMaxDiffP:
                findContact = False
                for k in range(length):
                    if k == i or k == j:
                        continue
                    if overlapA(tmpList[i], tmpList[k]) >= 0.9 and \
                            overlapA(tmpList[j], tmpList[k]) >= 0.9:
                        findContact = True
                        break
                if not findContact:
                    diffNum = diffNum + 1
            if overlapMax(tmpList[i], tmpList[j]) >= overMaxSameP and \
                    overlapMin(tmpList[i], tmpList[j]) <= overMinSameP and \
                    overlapA(tmpList[i], tmpList[j]) >= overMaxSameP:
                maxSameNum = maxSameNum + 1
        #判断第i个tanscript有没有可能成为分裂纠正的证据transcript,要求他尽可能完整,尽可能长,并且没有第三方连接证据
        if diffNum > 0 and maxSameNum == 0:
            if 'PS' in gene and len(gene) == 11:
                splitF = True
                break
    if not splitF:
        continue
    gtus = ttUFS(length+1)
    mgtu2MaxLen = [ tmpList[i].ed - tmpList[i].st + 1 for i in range(length) ]
    mgtu2MaxId = [i for i in range(length)]
    mgtu2TransNum=[getTranscriptFlncNum(tmpList[i].trans) for i in range(length)]
    for i in range(length):
        if tmpList[i].ed - tmpList[i].st + 1 >= maxTransSpan:
            continue
        for j in range(i + 1, length):
            if tmpList[j].ed - tmpList[j].st + 1 >= maxTransSpan:
                continue
            ir=gtus.find_set(i+1)-1
            jr=gtus.find_set(j+1)-1
            if overlapMax(tmpList[mgtu2MaxId[ir]],tmpList[mgtu2MaxId[jr]]) >=0.5:
                gtus.union(i+1, j+1)
                nid=0
                if mgtu2MaxLen[ir]>mgtu2MaxLen[jr]:
                    nid=mgtu2MaxId[ir]
                else:
                    nid=mgtu2MaxId[jr]
                nlen=max(mgtu2MaxLen[ir],mgtu2MaxLen[jr])
                mgtu2MaxLen[ir]=mgtu2MaxLen[jr]=nlen
                mgtu2MaxId[ir]=mgtu2MaxId[jr]=nid
                tmpTransNum=mgtu2TransNum[ir]+mgtu2TransNum[jr]
                mgtu2TransNum[ir]=mgtu2TransNum[jr]=tmpTransNum
    gtu2len = defaultdict(int)
    gtu2ind = defaultdict(int)
    for i in range(length):
        if tmpList[i].ed - tmpList[i].st + 1 >= maxTransSpan:
            continue
        gtu = gtus.find_set(i+1)
        tmpLen = tmpList[i].ed - tmpList[i].st + 1
        if gtu2len[gtu] < tmpLen:
            gtu2len[gtu] = tmpLen
            gtu2ind[gtu] = i
    splitGeneSet[gene] = 1
    for gtu in gtu2ind:
        if mgtu2TransNum[gtu-1]>1:
            spg2trans[gene].append(tmpList[gtu2ind[gtu]].trans)

# with open('data/spgPlot.dat', 'w') as of:
#     for gene in splitGeneSet:
#         if len(spg2trans[gene])<2:
#             continue
#         for trans in spg2trans[gene]:
#             print(gene, trans, sep='\t', file=of)

with open('data/splitGene.list.txt', 'w') as of:
    for gene in splitGeneSet:
        if len(spg2trans[gene])<2:
            continue
        print(gene, end='\t', file=of)
        for trans in spg2trans[gene]:
            print(trans, end='\t', file=of)
        print('', file=of)
# with open('data/spgGene.list', 'w') as of:
#     for gene in splitGeneSet:
#         if len(spg2trans[gene])<2:
#             continue
#         print(gene, file=of)
#输出供断点分析
for gene in gene2transBainfo:
    for i in gene2transBainfo[gene]:
        trans2Bainfo[i.trans]=i
with open('data/detect2breakSta.dat','w') as of:
    for gene in splitGeneSet:
        if len(spg2trans[gene])==2 and overlapMax(trans2Bainfo[spg2trans[gene][0]],trans2Bainfo[spg2trans[gene][1]])==0:
            bainfo1=trans2Bainfo[spg2trans[gene][0]]
            bainfo2=trans2Bainfo[spg2trans[gene][1]]
            print(gene,spg2trans[gene][0],spg2trans[gene][1],bainfo1.chr,min(bainfo1.ed,bainfo2.ed),max(bainfo1.st,bainfo2.st),sep='\t',file=of)
