import sys
sys.path.append('/home/xutun/transAs/src')
from ttlib.basicInfo import bainfo
from ttlib.basicInfo import calDis
import ttlib.geneStructure.nameServer as ns
import pandas as pd
class ep:
    def __init__(self,line):
        l = line.split('\t')
        self.ga = l[0]
        self.ea = int(l[1])
        self.pa = bainfo(l[2])
        self.gb = l[3]
        self.eb = int(l[4])
        self.pb = bainfo(l[5])

    def toStr(self):
        return f'{self.ga}\t{self.ea}\t{self.gb}\t{self.eb}'

    def simplePrint(self):
        print(self.ga,self.ea,self.gb,self.eb,sep='\t')


class epPair:
    def __init__(self,ep1,ep2):
        self.ep1 = ep1
        self.ep2 = ep2
        self.disa = calDis(ep1.pa,ep2.pa)
        self.disb = calDis(ep1.pb,ep2.pb)
        self.diffDis = self.disa - self.disb
        self.ena = self.ep1.ea - self.ep2.ea
        self.enb = self.ep1.eb - self.ep2.eb
        if self.ena < 0:
            self.ena = -self.ena
            self.enb = -self.enb
        self.diffEn = self.ena - self.enb
    def toStr(self):
        return f'{self.ep1.ga}\t{self.ep1.gb}\t{self.ena}\t{self.enb}\t{self.diffEn}\t{self.disa}\t{self.disb}\t{self.diffDis}\t{self.ep1.ea}\t{self.ep2.ea}\t{self.ep1.eb}\t{self.ep2.eb}'

class epChain:
    def __init__(self,epList):
        self.eps = epList

    def geneAdjEpPair(self):
        relList = []
        for i in range(len(self.eps)-1):
            relList.append(epPair(self.eps[i],self.eps[i+1]))
        return relList

    def toStr(self):
        retStr = ''
        for ep in self.eps:
            retStr = retStr + ep.toStr() + '\n'
        return retStr

    def print(self):
        for ep in self.eps:
            ep.simplePrint()

