import os
import sys
sys.path.append('/home/xutun/transAs')
from ttlib.stringBase import  generateWorkD
import ttlib.geneStructure.nameServer as ns
from ttlib.geneStructure.dataStruct import ep,epChain
from ttlib.geneStructure.relaServer import getMaxIncreaseSeq,getInvSeq,getDupSeq

epChainList = []
invChainList = []
dupChainList = []

def geneEpChain(tupleList,epList):
    tupleList.sort()
    tmpList = []
    for tup in tupleList:
        for ep in epList:
            if ep.ea == tup[0] and ep.eb == tup[1]:
                tmpList.append(ep)
                break
    return epChain(tmpList)

def getExonTupleList(epList):
    relList = []
    for ep in epList:
        relList.append((ep.ea,ep.eb))
    return relList

def solveCache(cache,_epChainList):
    global invChainList
    if len(cache) == 0:
        return
    epList = []
    for line in cache:
        epList.append(ep(line))
    if epList[0].ga == epList[0].gb:
        return
    epInfo = getExonTupleList(epList)

    #normalDetect
    maxNum,tupleList = getMaxIncreaseSeq(epInfo)
    if maxNum > 1:
        _epChainList.append(geneEpChain(tupleList,epList))

    #invDetect
    # maxNum, tupleList = getInvSeq(epInfo)
    # if maxNum > 2:
    #     tmpChain = geneEpChain(tupleList, epList)
    #     invChainList.append(tmpChain)

    #dupDetect
    maxNum,tupleList = getDupSeq(epInfo)
    if maxNum > 2:
        tmpChain = geneEpChain(tupleList,epList)
        dupChainList.append(tmpChain)

def geneStructDetect(filteredBlast,workDir=''):
    if workDir == '':
        workDir = generateWorkD()
    os.system(f'mkdir -p {workDir}')
    pairInfoF = f'{workDir}/pairInfo.tab'
    ns.reformatBlast(filteredBlast,pairInfoF)
    prePair = ''
    cache = []

    for line in open(pairInfoF,'r'):
        line = line.strip()
        l = line.split('\t')
        nowPair = l[0] + l[3]
        if nowPair != prePair:
            solveCache(cache,epChainList)
            cache = []
            prePair  = nowPair
        cache.append(line)
    solveCache(cache, epChainList)

    with open(f'{workDir}/dupChain.tab','w') as of:
        for ch in dupChainList:
            print(ch.toStr(),file=of,end='')

    # with open(f'{workDir}/invChain.tab','w') as of:
    #     for ch in invChainList:
    #         print(ch.toStr(),file=of,end='')

    with open(f'{workDir}/normalChain.tab','w') as of:
        for ch in epChainList:
            print(ch.toStr(),file=of,end='')

    with open(f'{workDir}/epPair.tab','w') as of:
        print('GA\tGB\tENA\tENB\tDIFFEN\tDISA\tDISB\tDIFFDIS\tGAEA\tGAEB\tGBEA\tGBEB',file=of)
        for epChain in epChainList:
            for epPair in epChain.geneAdjEpPair():
                print(epPair.toStr(),file = of)
