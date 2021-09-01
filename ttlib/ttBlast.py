import os
from collections import defaultdict
import sys
sys.path.append('/data/ttlib')
sys.path.append('/home/xutun/transAs/src/ttlib')
sys.path.append('/data/home/xutun/mySrc/ttlib')
from readBase import getId2lenFromFai
from ttWriter import ttWriter
# samtoolsP = '/home/xutun/miniconda3/envs/cppEnv/bin/samtools'
samtoolsP = '/data/home/xutun/miniconda3/envs/tt/bin/samtools'
def blast(toolName,faa,fab,workd,extraPara=''):
    # toolD='/home/xutun/miniconda3/envs/norm/bin'
    # toolD='/home/xutun/miniconda3/envs/cppEnv/bin'
    toolD='/data/home/xutun/miniconda3/envs/tt/bin'
    faaName=faa.split('/')[-1]
    fabName=fab.split('/')[-1]
    dbtype='prot'
    if toolName=='blastn':
        dbtype='nucl'
    os.system("""
    if [ -d %s ];then
        rm -rf %s
    fi
    mkdir %s
    """%(workd,workd,workd))
    os.system("""
    %s/makeblastdb -in %s -dbtype %s -parse_seqids -out %s/%s.db
    """%(toolD,faa,dbtype,workd,faaName))
    os.system("""
    %s/%s -query %s -db %s/%s.db -out %s/%s.%s.out -evalue 1e-3 -num_threads 48 -outfmt 6 -num_alignments 100 %s
    """%(toolD,toolName,fab,workd,faaName,workd,faaName,fabName,extraPara))
    return """%s/%s.%s.out"""%(workd,faaName,fabName)

def filterBlastRelCheckCache(cacheList,aId2len,bId2len,overlapA,overlapB,identity,writer):
    if len(cacheList) == 0:
        return
    aSeq = []
    bSeq = []
    aId = ''
    bId = ''
    for line in cacheList:
        line = line.split('\t')
        aId = line[1]
        bId = line[0]
        if float(line[2]) >= identity:
            aSeq.append(( min(int(line[8]),int(line[9])),max(int(line[8]),int(line[9])) ))
            bSeq.append(( min(int(line[6]),int(line[7])),max(int(line[6]),int(line[7])) ))
    #merge adjacent seq
    totLen = len(aSeq)
    i = 0
    findMerge = False
    while i < totLen:
        j = i + 1
        while j < totLen:
            if ( aSeq[i][0] >= aSeq[j][0] and aSeq[i][0] <= aSeq[j][1] ) or \
               ( aSeq[i][1] >= aSeq[j][0] and aSeq[i][1] <= aSeq[j][1]):
                aSeq[i] = ( min(aSeq[i][0],aSeq[j][0]),max(aSeq[i][1],aSeq[j][1]) )
                aSeq[j] = aSeq[-1]
                aSeq.pop()
                findMerge = True
                totLen = totLen - 1
                break
            else:
                j = j + 1
        if findMerge:
            i = 0
            findMerge = False
        else:
            i= i + 1
    totLen = len(bSeq)
    i = 0
    findMerge = False
    while i < totLen:
        j = i + 1
        while j < totLen:
            if (bSeq[i][0] >= bSeq[j][0] and bSeq[i][0] <= bSeq[j][1]) or \
                    (bSeq[i][1] >= bSeq[j][0] and bSeq[i][1] <= bSeq[j][1]):
                bSeq[i] = (min(bSeq[i][0], bSeq[j][0]), max(bSeq[i][1], bSeq[j][1]))
                bSeq[j] = bSeq[-1]
                bSeq.pop()
                findMerge = True
                totLen = totLen - 1
                break
            else:
                j = j + 1
        if findMerge:
            i = 0
            findMerge = False
        else:
            i = i + 1
    aRange = 0
    bRange = 0
    for seq in aSeq:
        aRange = aRange + (seq[1]-seq[0]+1)
    for seq in bSeq:
        bRange = bRange + (seq[1]-seq[0]+1)
    if aRange/aId2len[aId] >= overlapA and bRange/bId2len[bId] >= overlapB:
        for line in cacheList:
            cpline = line
            line = line.split('\t')
            if float(line[2]) >= identity:
                writer.writ(cpline)

def filterBlastRel(resultF,fileA,fileB,overlapA,overlapB,identity,outF):
    myWriter = ttWriter(outF)
    os.system(f'{samtoolsP} faidx {fileA} -o {fileA}.fai')
    os.system(f'{samtoolsP} faidx {fileB} -o {fileB}.fai')
    aId2len = getId2lenFromFai(f'{fileA}.fai')
    print(f'{fileA}.fai')
    bId2len = getId2lenFromFai(f'{fileB}.fai')
    os.system(f'sort {resultF} -k1,1 -k2,2 >{resultF}.sort')
    cacheList = []
    nowTotId = ''
    with open(f'{resultF}.sort','r') as f:
        line=f.readline().strip()
        while line != '':
            cpline = line
            line = line.split('\t')
            aId = line[1]
            bId = line[0]
            totId = aId + bId
            if totId != nowTotId:
                filterBlastRelCheckCache(cacheList,aId2len,bId2len,overlapA,overlapB,identity,myWriter)
                cacheList = []
                nowTotId = totId
            cacheList.append(cpline)
            line = f.readline().strip()
    myWriter.close()

def filterReverseBackAlign(resultF,fileA,fileB,overlapA,overlapB,identity,outF):
    myWriter = ttWriter(outF)
    os.system(f'{samtoolsP} faidx {fileA} -o {fileA}.fai')
    os.system(f'{samtoolsP} faidx {fileB} -o {fileB}.fai')
    aId2len = getId2lenFromFai(f'{fileA}.fai')
    bId2len = getId2lenFromFai(f'{fileB}.fai')
    # os.system(f'sort {resultF} -k1,1 -k2,2 >{resultF}.sort')
    cacheList = []
    nowTotId = ''
    tmpNum = 0
    with open(f'{resultF}', 'r') as f:
        line = f.readline().strip()
        while line != '':
            cpline = line
            line = line.split('\t')
            aId = line[1]
            bId = line[0]
            dis1 = int(line[7]) - int(line[6])
            dis2 = int(line[9]) - int(line[8])
            if ((dis1<0 and dis2>0) or (dis1>0 and dis2<0)) and aId==bId:
                totId = aId + bId
                if totId != nowTotId:
                    filterBlastRelCheckCache(cacheList, aId2len, bId2len, overlapA, overlapB, identity, myWriter)
                    cacheList = []
                    nowTotId = totId
                    tmpNum = tmpNum + 1
                cacheList.append(cpline)
            line = f.readline().strip()
    print(f'tmpNum:{tmpNum}')
    myWriter.close()