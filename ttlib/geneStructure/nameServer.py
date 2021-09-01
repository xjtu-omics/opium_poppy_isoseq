import sys
import os
import re
sys.path.append('/data/ttlib')
from ttlib.basicInfo import bainfo,calDis

#To blast format
#DCW010001.0++Exon2::chr1:12753-13489(+)
#DCW010001.0++Intron1::chr1:12524-12753(+)
#{geneName}++[Exon|Intron]{Number}::{seqName}:{st}-{ed}({strand})

def getGeneName(blastId):
    return blastId.split('++')[0]

def getExonName(blastId):
    return blastId.split('::')[0]

def getExonNum(blastId):
    return int(re.findall('\+\+.*([0-9]+)::',blastId)[0])

def getBasicInfoStr(blastId):
    return blastId.split('::')[1]

def getBasicInfo(blastId):
    chr = re.findall('::(.*)?:',blastId)[0]
    st = int(re.findall(':([0-9]+)?\-',blastId)[0])
    ed = int(re.findall('\-([0-9]+)?\(',blastId)[0])
    strand = re.findall('\((.{1})\)',blastId)[0]
    return bainfo(chr,st,ed,blastId,strand)

def getRefDisFromBl(ida,sta,eda,idb,stb,edb):
    apos = getBasicInfo(ida)
    bpos = getBasicInfo(idb)
    if apos.strand == '+':
        apos.ed = apos.st + eda
        apos.st = apos.st + sta
    else:
        apos.st = apos.ed - eda
        apos.ed = apos.ed - sta
    if bpos.strand == '+':
        bpos.ed = bpos.st + edb
        bpos.st = bpos.st + stb
    else:
        bpos.st = bpos.ed - edb
        bpos.ed = bpos.ed - stb
    return calDis(apos,bpos)

def getRefDis(bla,blb):
    bla = bla.split('\t')
    blb = blb.split('\t')
    refEDis = getRefDisFromBl(bla[0],int(bla[6]),int(bla[7]),
                              blb[0],int(blb[6]),int(blb[7]))
    tarEDis = getRefDisFromBl(bla[1],int(bla[8]),int(bla[9]),
                              blb[1],int(blb[8]),int(blb[9]))
    return bla[0],bla[1],refEDis,tarEDis

def getExonRefDis(bla,blb):
    bla = bla.split('\t')
    blb = blb.split('\t')
    refPos0 = getBasicInfo(bla[0])
    refPos1 = getBasicInfo(blb[0])
    tarPos0 = getBasicInfo(bla[1])
    tarPos1 = getBasicInfo(blb[1])
    return bla[0],bla[1],blb[0],blb[1],calDis(refPos0,refPos1),calDis(tarPos0,tarPos1)

def reformatBlast(blastf,pairf):
    with open(pairf,'w') as of:
        for line in open(blastf,'r'):
            line = line.split('\t')
            geneA = getGeneName(line[0])
            geneB = getGeneName(line[1])
            exonA = getExonNum(line[0])
            exonB = getExonNum(line[1])
            ginfoA = getBasicInfoStr(line[0])
            ginfoB = getBasicInfoStr(line[1])
            print(f'{geneA}\t{exonA}\t{ginfoA}\t{geneB}\t{exonB}\t{ginfoB}',file=of)
    os.system(f'sort -k1,1 -k4,4 {pairf} >tmp && mv tmp {pairf}')