import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import math
import sys
sys.path.append('/data/mySrc')
sys.path.append('/data/home/xutun/mySrc')
from ttlib.ttplot import gepardDotplot
from ttlib.basicInfo import ttSeq
class ttdp:
    def __init__(self,raname,raseq='',rbname='',rbseq=''):
        if raseq != '' and rbname != '' and rbseq != '':
            self.an = raname
            self.bn = rbname
            self.aseq = str(raseq)
            self.bseq = str(rbseq)
            self.alen = len(raseq)
            self.blen = len(rbseq)
            self.dp = [[0 for i in range(self.blen)] for j in range(self.alen)]
        if raseq == '' and rbname == '' and rbseq == '':
            seqTup = raname
            self.an = seqTup[0].qn
            self.bn = seqTup[1].qn
            self.aseq = seqTup[0].seq
            self.bseq = seqTup[1].seq
            self.alen = len(self.aseq)
            self.blen = len(self.bseq)
            self.dp = [[0 for i in range(self.blen)] for j in range(self.alen)]

    def generateDp(self,bin=10,cutoff=9,strand='+'):
        self.dp = [[0 for i in range(self.blen)] for j in range(self.alen)]
        for ast in range(bin-1,self.alen):
            sameNum = 0
            for i in range(bin-1):
                sameNum = sameNum + (self.aseq[i+ast-bin+1] == self.bseq[i])
            i = ast
            j = bin-1
            while i < self.alen and j < self.blen:
                if self.aseq[i] == self.bseq[j]:
                    sameNum = sameNum + 1
                if sameNum >= cutoff:
                    self.dp[i][j] = 1
                if self.aseq[i-bin+1] == self.bseq[j-bin+1]:
                    sameNum = sameNum - 1
                i = i+1
                j = j+1
        for bst in range(bin - 1, self.blen):
            sameNum = 0
            for i in range(bin - 1):
                sameNum = sameNum + (self.bseq[i + bst - bin + 1] == self.aseq[i])
            i = bin - 1
            j = bst
            while i < self.alen and j < self.blen:
                if self.aseq[i] == self.bseq[j]:
                    sameNum = sameNum + 1
                if sameNum >= cutoff:
                    self.dp[i][j] = 1
                if self.aseq[i - bin + 1] == self.bseq[j - bin + 1]:
                    sameNum = sameNum - 1
                i = i + 1
                j = j + 1
        return self.dp

    def showHoughTrans(self):
        tots = []
        totr = []
        for i in range(self.alen):
            for j in range(self.blen):
                if self.dp[i][j] == 1:
                    for theta in range(130,141):
                        s = theta * math.pi / 180
                        r = i * math.cos(s) + j * math.sin(s)
                        tots.append(s)
                        totr.append(r)
        sns.jointplot(x=tots,y=totr,kind='kde')
        plt.show()
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.scatter(x=tots,y=totr,s=0.1)
        plt.show()
        plt.close()

    def plotDp(self,fileName='tmp.png'):
        #gepardDotplot(self.an,self.aseq,self.bn,self.bseq,'tmp.png')
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        x = []
        y = []
        for i in range(self.alen):
            for j in range(self.blen):
                if self.dp[i][j] == 1:
                    x.append(i)
                    y.append(j)
        print(len(x))
        ax.scatter(x,y,s=0.1,alpha=0.1)
        plt.savefig(fileName)
