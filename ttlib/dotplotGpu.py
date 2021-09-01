import pandas as pd
import numpy as np
import pycuda.autoinit
import pycuda.driver as drv
from pycuda import gpuarray
from pycuda.compiler import SourceModule
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import collections as mc
from matplotlib import colors as mcolors
import math
from operator import attrgetter
import sys
from timeit import default_timer as timer
sys.path.append('/data/mySrc')
sys.path.append('/data/home/xutun/mySrc')
sys.path.append('/public/xutun2/mySrc')
from ttlib.ttplot import gepardDotplot
from ttlib.basicInfo import ttSeq
from ttlib.gpuBase import getCKerFunction
from ttlib.gpuBase import sumGpu
from ttlib.ttDataStructure import ttUFS
kerP = '/data/home/xutun/mySrc/ttlib/ker'
#kerP = '/public/xutun2/mySrc/ttlib/ker'

class alignSegment:
    def __init__(self,st,ed,id):
        self.st = st
        self.ed = ed
        self.id = id
        self.intercept = self.st[0] - self.st[1]

class ttdpGpu:
    def __init__(self,raname,myId,raseq='',rbname='',rbseq=''):
        # !!!To change
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

            self.aseq = []
            for nucl in seqTup[0].seq:
                self.aseq.append(ord(nucl))
            self.aseq = np.array(self.aseq,dtype=np.int8)

            self.bseq = []
            for nucl in seqTup[1].seq:
                self.bseq.append(ord(nucl))
            self.bseq = np.array(self.bseq, dtype=np.int8)

            self.alen = np.int32(min(5000,len(self.aseq)))
            self.blen = np.int32(min(5000,len(self.bseq)))
            # self.alen = np.int32(len(self.aseq))
            # self.blen = np.int32(len(self.bseq))
            self.dp = np.zeros(self.alen * self.blen,dtype=np.int8)
            self.relGpu = None
            self.myId = myId

    def generateHeatDp(self,minKmer=4,maxKmer=33):
        heatDpNvcc = getCKerFunction(f'{kerP}/dotPlotSvKer.cpp','generateDp')
        N = np.int32(self.alen*self.blen)
        heatdp = np.zeros(N,dtype=np.int8)
        for i in range(minKmer,maxKmer+1,1):
            heatDpNvcc(drv.In(self.aseq),drv.In(self.bseq),drv.InOut(heatdp),self.alen,self.blen,N,np.int32(i),np.int32(i),block=(512,1,1),grid=(int(np.ceil(N/512)),1))
        heatdp = np.array(heatdp).reshape(self.blen,self.alen)
        fig,ax=plt.subplots()
        heatdp[heatdp==0] = 1
        sns.heatmap(data=np.log2(heatdp),cbar=False,xticklabels=False,yticklabels=False)
        fig.set_size_inches(40,40)
        fig.tight_layout(pad=0)
        plt.savefig('plot/heatPlot.png')
        plt.close()
        exit()
        subSize=50
        sbpx = int(self.alen/subSize)
        sbpy = int(self.blen/subSize)
        fig,axes = plt.subplots(sbpx,sbpy,figsize=(40,40))
        for i in range(sbpx):
            for j in range(sbpy):
                tpd = heatdp[i*subSize:i*subSize+subSize,j*subSize:j*subSize+subSize]
                print(i,j)
                sns.heatmap(ax=axes[i,j],data=tpd,cbar=False,xticklabels=False,yticklabels=False)
        fig.tight_layout(pad=0)
        plt.savefig('plot/multiHeatPlot.png',pad_inches=0,bbox_inches='tight')

    def _1to2(self,i):
        return [i%self.alen,int(i/self.alen)]

    def _2to1(self,pos):
        return pos[1]*self.alen + pos[0]

    def testShowNowSegments(self,N,outP):
        st = self.stPos.get()
        ed = self.edPos.get()
        oneDimensionLines = []
        blastLines = []
        segNum = 0
        for i in range(N):
            if st[i]==i:
                segNum = segNum + 1
                oneDimensionLines.append([i,ed[i]])
        for odl in oneDimensionLines:
            blastLines.append([self._1to2(odl[0]),self._1to2(odl[1])])
        fig, ax = plt.subplots()
        lc = mc.LineCollection(blastLines)
        ax.add_collection(lc)
        ax.autoscale()
        fig.set_size_inches(40*(self.alen/max(self.alen,self.blen)), 40*(self.blen/max(self.alen,self.blen)))
        plt.savefig(outP)
        plt.close()

    def generateProbDp(self,maxBin):
        drv.init()
        dev = drv.Device(2)
        ctx = dev.make_context()
        aseqGpu = gpuarray.to_gpu(self.aseq)
        bseqGpu = gpuarray.to_gpu(self.bseq)
        # print(self.alen,self.blen)
        # print('tot',self.alen*self.blen)
        N = np.int32(self.alen*self.blen)
        self.probRelGpu = gpuarray.to_gpu(np.zeros(N,dtype=np.float32))
        self.stPos = gpuarray.to_gpu(np.zeros(N,dtype=np.int32))
        self.edPos = gpuarray.to_gpu(np.zeros(N,dtype=np.int32))
        tmpStPos = gpuarray.to_gpu(np.zeros(N,dtype=np.int32))
        tmpEdPos = gpuarray.to_gpu(np.zeros(N,dtype=np.int32))
        tmpFlag = gpuarray.to_gpu(np.zeros(N,dtype=np.int32))
        kerGnerateProbDp = getCKerFunction(f'{kerP}/dotPlotSvKer.cpp','generateProbDp')
        kerInitStEd = getCKerFunction(f'{kerP}/dotPlotSvKer.cpp','initStEd')
        kerGetForwardRela = getCKerFunction(f'{kerP}/dotPlotSvKer.cpp','findPairRelaForward')
        kerGetBackwardRela = getCKerFunction(f'{kerP}/dotPlotSvKer.cpp','findPairRelaBackward')
        kerFindPairPair = getCKerFunction(f'{kerP}/dotPlotSvKer.cpp','findPairPair')
        kerUpdateStEd = getCKerFunction(f'{kerP}/dotPlotSvKer.cpp','updateStEd')
        kerStaSegNum = getCKerFunction(f'{kerP}/dotPlotSvKer.cpp','staSegNum')
        kerInitArray = getCKerFunction(f'{kerP}/gpuBaseKer.cpp','initArray')
        kerGnerateProbDp(aseqGpu,bseqGpu,self.probRelGpu,self.alen,self.blen,N,np.int32(maxBin),block=(1024,1,1),grid=(int(np.ceil(N/1024)),1))
        # self.plotHeatDotplot()
        kerInitStEd(self.probRelGpu,self.stPos,self.edPos,N,block=(1024,1,1),grid=(int(np.ceil(N/1024)),1))
        self.initStFlag()
        while True:
            st = timer()
            kerInitArray(tmpStPos,np.int32(-1),N,block=(1024,1,1),grid=(int(np.ceil(N/1024)),1))
            kerInitArray(tmpEdPos,np.int32(-1),N,block=(1024,1,1),grid=(int(np.ceil(N/1024)),1))
            kerInitArray(tmpFlag,np.int32(0),N,block=(1024,1,1),grid=(int(np.ceil(N/1024)),1))

            # kerInitArray(tmpFlag,np.int32(0),N,block=(1024,1,1),grid=(int(np.ceil(N/1024)),1))
            kerGetForwardRela(tmpFlag,self.segStPos,self.stPos,self.edPos,tmpEdPos,self.alen,self.blen,self.segNum,block=(1024,1,1),grid=(int(np.ceil(self.segNum/1024)),1))
            # kerGetForwardRela(tmpFlag,self.segStPos,self.stPos,self.edPos,tmpEdPos,self.alen,self.blen,N,block=(1024,1,1),grid=(int(np.ceil(N/1024)),1))
            # print('tmpFlag:::',np.sum(tmpEdPos.get()>=0))
            # exit()

            # kerInitArray(tmpFlag,np.int32(0),N,block=(1024,1,1),grid=(int(np.ceil(N/1024)),1))
            # kerGetBackwardRela(tmpFlag,self.segStPos,self.stPos,self.edPos,tmpStPos,self.alen,self.blen,N,block=(1024,1,1),grid=(int(np.ceil(N/1024)),1))
            kerGetBackwardRela(tmpFlag,self.segStPos,self.stPos,self.edPos,tmpStPos,self.alen,self.blen,self.segNum,block=(1024,1,1),grid=(int(np.ceil(self.segNum/1024)),1))
            # print('tmpFlag:::',np.sum(tmpFlag.get()))
            # exit()


            # kerGetBackwardRela(self.segStPos,self.stPos,self.edPos,tmpStPos,self.alen,self.blen,self.segNum,block=(1024,1,1),grid=(int(np.ceil(self.segNum/1024)),1))
            # kerGetForwardRela(self.segStPos,self.stPos,self.edPos,tmpEdPos,self.alen,self.blen,N,block=(1024,1,1),grid=(int(np.ceil(N/1024)),1))
            # kerGetBackwardRela(self.segStPos,self.stPos,self.edPos,tmpStPos,self.alen,self.blen,N,block=(1024,1,1),grid=(int(np.ceil(N/1024)),1))
            # print('c')
            kerFindPairPair(tmpStPos,tmpEdPos,tmpFlag,self.stPos,self.edPos,N,block=(1024,1,1),grid=(int(np.ceil(N/1024)),1))
            newRelaNum = sumGpu(tmpFlag,N)
            # print('d')
            kerUpdateStEd(self.stPos,self.edPos,self.segStFlag,self.segStPos,self.segNum,N,block=(1024,1,1),grid=(int(np.ceil(N/1024)),1))
            # print('f')
            kerStaSegNum(self.stPos,tmpFlag,N,block=(1024,1,1),grid=(int(np.ceil(N/1024)),1))
            # print('g')
            segNum = sumGpu(tmpFlag,N)
            # print('flagSum:',np.sum(self.segStFlag.get()))
            # print(f'GPU Cal:{timer()-st}\tsegNumj:{segNum}\tnewRelaNum:{newRelaNum}')
            # print('self.segNum:',self.segNum)
            self.updateSeg()
            if newRelaNum==0:
                # self.testShowNowSegments(N,f'plot/{newRelaNum}.png')
                break
        #Fix the end of the segment
        kerFixCheck = getCKerFunction(f'{kerP}/dotPlotSvKer.cpp','fixCheck')
        self.initFixData()
        while len(self.stSearchPosCpu)>0:
            # print('To Search num:',len(self.stSearchPosCpu))
            kerFixCheck(self.stSearchPos,self.probRelGpu,self.fixRel,np.int32(self.alen*self.blen),self.alen,self.blen,np.int32(len(self.stSearchPosCpu)),block=(1024,1,1),grid=(int(np.ceil(N/1024)),1))
            self.updateFixData()
        self.mergeOverlapReads()

        #记录断点修复后的结果方便作图
        self._2StepSegStPosCpu = np.array(self.segStPosCpu,dtype=np.int32)
        self._2StepSegEdPosCpu = np.array(self.segEdPosCpu,dtype=np.int32)

        self.connectSegments()

        self._3StepSegStPosCpu = self.segStPosCpu
        self._3StepSegEdPosCpu = self.segEdPosCpu

        nowTime = timer()

        self.plotSegment()
        ctx.pop()
        return nowTime

        minSpecity = 5


        return
        fig, ax = plt.subplots()
        # blastLines = [[(886,724),(1966,1800)],[(1,14),(498,511)],[(554,503),(787,734)]]
        # blastLines = [[(198,887),(1114,1839)],[(507,1395),(692,1590)],[(721,1237),(865,1382)]]
        # blastLines = [[(189,134),(918,865)],[(1,12),(588,591)],[(130,339),(404,617)]]
        blastLines = [[(1533,1447),(3244,3158)],[(5,1),(1447,1446)],[(1308,379),(2106,1155)],
                      [(610,1465),(1156,2020)],[(1187,174),(1233,220)],[(175,1186),(221,1232)]]
        lc = mc.LineCollection(blastLines,colors=np.array([(1, 0, 0, 1), (0, 1, 0, 1), (0, 0, 1, 1)]))
        # ax.add_collection(lc)
        probdp[probdp<=1e-6] = 0
        probdp[probdp>=1e-6] = 1
        probdp = 1-probdp
        sns.heatmap(data=probdp, cbar=False, xticklabels=False, yticklabels=False)
        fig.set_size_inches(40*(self.alen/max(self.alen,self.blen)), 40*(self.blen/max(self.alen,self.blen)))
        fig.tight_layout(pad=0)
        plt.savefig('plot/probPlot.png')
        plt.close()
        sns.kdeplot(-np.log10(np.array(probdp).reshape(self.alen*self.blen,)),shade=True)
        fig.set_size_inches(40,40)
        plt.savefig('plot/kdePlot.png')
        plt.close()
        ctx.pop()

    def initFixData(self):
        self.fcCutOff = 1/1.05
        self.checkLen = 5
        self.checkIt = 5
        self.checkStep = 3
        self.segStPosCpu = self.segStPos.get()
        self.edPosCpu = self.edPos.get()
        self.segEdPosCpu = np.array(self.edPosCpu[self.segStPosCpu],dtype=np.int32)
        self.segNum = len(self.segEdPosCpu)
        #With dynamic length since some endPos will be asigned inactive
        self.stSearchPosCpu = np.zeros(len(self.segStPosCpu)*2*self.checkIt,dtype=np.int32)

        #为作图做准备
        self._1StepSegStPosCpu = np.array(self.segStPosCpu.tolist(),dtype=np.int32)
        self._1StepSegEdPosCpu = np.array(self.segEdPosCpu.tolist(),dtype=np.int32)


        stSearchIdx = 0
        for st in self.segStPosCpu:
            for adjustPos in range(self.checkIt):
                _2Pos = self._1to2(st)
                _2Pos[0] = _2Pos[0]-1-self.checkStep*adjustPos
                _2Pos[1] = _2Pos[1]-1-self.checkStep*adjustPos
                toSearch1Pos = self._2to1(_2Pos)
                self.stSearchPosCpu[stSearchIdx] = toSearch1Pos
                stSearchIdx = stSearchIdx+1
        for ed in self.segEdPosCpu:
            for adjustPos in range(self.checkIt):
                _2Pos = self._1to2(ed)
                _2Pos[0] = _2Pos[0]+self.checkLen+self.checkStep*adjustPos
                _2Pos[1] = _2Pos[1]+self.checkLen+self.checkStep*adjustPos
                toSearch1Pos = self._2to1(_2Pos)
                self.stSearchPosCpu[stSearchIdx] = toSearch1Pos
                stSearchIdx = stSearchIdx+1

        #flagSt1,flagSt2, ... ,flagStN,flagEd1,flagEd2, ... ,flagEdN
        self.activeEndPosCpu = np.ones(len(self.segStPosCpu)*2)

        self.stSearchPos = gpuarray.to_gpu(self.stSearchPosCpu)
        self.fixRelCpu = np.zeros(len(self.stSearchPosCpu),dtype=np.float32)
        self.fixRel = gpuarray.to_gpu(self.fixRelCpu)


        #为作图做准备

    def plotFixCheckRel(self):
        self.fixRelCpu = self.fixRel.get()
        self.fixRelCpu = np.array(self.fixRelCpu[self.fixRelCpu!=1])
        pltd = pd.DataFrame({'fixRel':self.fixRelCpu})
        sns.displot(pltd,x='fixRel',kind='kde')
        plt.savefig('plot/fixRelDensity.pdf')
        plt.close()

    def initStFlag(self):
        tmpSt = self.stPos.get()
        self.segStFlag = np.zeros(len(tmpSt),dtype=np.int32)
        self.segStFlag[tmpSt==np.arange(len(tmpSt))] = 1
        self.segStPos = np.array(np.arange(len(tmpSt))[np.array(self.segStFlag,dtype=bool)],dtype=np.int32)
        self.segNum = np.int32(len(self.segStPos))
        self.segStFlag = gpuarray.to_gpu(np.zeros(len(self.segStPos),dtype=np.int32))
        self.segStPos = gpuarray.to_gpu(self.segStPos)

    def updateSeg(self):
        tmpStFlag = self.segStFlag.get()
        tmpSegStPos = self.segStPos.get()
        tmpSegStPos = np.array(tmpSegStPos[np.array(tmpStFlag,dtype=bool)],dtype=np.int32)
        self.segNum = np.int32(len(tmpSegStPos))
        self.segStPos = gpuarray.to_gpu(tmpSegStPos)
        self.segStFlag = gpuarray.to_gpu(np.zeros(len(tmpSegStPos),dtype=np.int32))

    def plotHeatDotplot(self):
        fig, ax = plt.subplots()
        probdp = self.probRelGpu.get()
        tmpMin = np.min(probdp[probdp > 0])
        probdp[probdp <= 0] = tmpMin / 10
        sns.heatmap(data=np.array(-np.log10(probdp)).reshape(self.blen, self.alen),
                    cbar=False, xticklabels=False,yticklabels=False)
        fig.set_size_inches(40 * (self.alen / max(self.alen, self.blen)), 40 * (self.blen / max(self.alen, self.blen)))
        fig.tight_layout(pad=0)
        plt.savefig('plot/probPlott.png')
        plt.close()

    def generateDp(self,bin=10,cutoff=9,strand='+'):
        kerCppf = open(f'{kerP}/generateDp.cpp','r')
        mod = SourceModule(kerCppf.read())
        kerCppf.close()
        cexec = mod.get_function('generateDp')
        N = np.int32(self.alen*self.blen)
        cexec(drv.In(self.aseq),drv.In(self.bseq),drv.Out(self.dp),self.alen,self.blen,N,np.int32(10),np.int32(10),block=(512,1,1),grid=(int(np.ceil(self.alen*self.blen/512)),1))
        fig, ax = plt.subplots()
        sns.heatmap(data=np.array(self.dp).reshape(self.blen,self.alen), cbar=False, xticklabels=False, yticklabels=False)
        fig.set_size_inches(40, 40)
        fig.tight_layout(pad=0)
        plt.savefig('plot/dotPlot.png')
        plt.close()
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
        print(sum(self.dp))
        dotIdx = np.where(self.dp)
        x = dotIdx / self.alen
        y = dotIdx % self.alen
        ax.scatter(x,y,s=0.1,alpha=0.1)
        plt.savefig(fileName)

    def updateFixData(self):
        self.fixRelCpu = self.fixRel.get()
        self.stSearchPosCpu = []
        posIdx = 0
        calRelIdx = 0
        tmpToSearchSt = []
        tmpToSearchEd = []
        #检查和更新所有候选点的迭代情况
        #一半大循环，对所有的线段起点进行检查
        while posIdx<len(self.segStPosCpu):
            if self.activeEndPosCpu[posIdx] == 0:
                posIdx = posIdx+1
                continue
            #以下是该位置在有在上次的迭代搜索的范围内的情况
            preFailNum = 0
            maxExtendLen = 0
            for i in range(self.checkIt):
                if preFailNum > 2:
                    calRelIdx = calRelIdx+1
                    continue
                #以下是在该断点有连续比对的情况
                if self.fixRelCpu[calRelIdx] > self.fcCutOff:
                    preFailNum = preFailNum + 1
                else:
                    preFailNum = 0
                if preFailNum == 0:
                    maxExtendLen = i+1
                calRelIdx = calRelIdx+1
            if maxExtendLen>0:
                _2Pos = self._1to2(self.segStPosCpu[posIdx])
                tmp2Pos = _2Pos
                _2Pos[0] = _2Pos[0] - maxExtendLen*self.checkStep
                _2Pos[1] = _2Pos[1] - maxExtendLen*self.checkStep
                _1Pos = self._2to1(_2Pos)
                self.segStPosCpu[posIdx] = _1Pos
                tmpToSearchSt.append(_1Pos)
            else:
                self.activeEndPosCpu[posIdx] = 0
            posIdx = posIdx+1
        #另一半大循环，对所有的线段终点进行检查
        while posIdx < 2*len(self.segStPosCpu):
            if self.activeEndPosCpu[posIdx] == 0:
                posIdx = posIdx + 1
                continue
            # 以下是该位置在有在上次的迭代搜索的范围内的情况
            preFailNum = 0
            maxExtendLen = 0
            for i in range(self.checkIt):
                if preFailNum > 2:
                    calRelIdx = calRelIdx + 1
                    continue
                # 以下是在该断点有连续比对的情况
                if self.fixRelCpu[calRelIdx] > self.fcCutOff:
                    preFailNum = preFailNum + 1
                else:
                    preFailNum = 0
                if preFailNum == 0:
                    maxExtendLen = i + 1
                calRelIdx = calRelIdx + 1
            if maxExtendLen > 0:
                _2Pos = self._1to2(self.segEdPosCpu[posIdx-len(self.segStPosCpu)])
                tmp2Pos = _2Pos
                _2Pos[0] = _2Pos[0] + maxExtendLen * self.checkStep
                _2Pos[1] = _2Pos[1] + maxExtendLen * self.checkStep
                _1Pos = self._2to1(_2Pos)
                self.segEdPosCpu[posIdx-len(self.segStPosCpu)] = _1Pos
                tmpToSearchEd.append(_1Pos)

            else:
                self.activeEndPosCpu[posIdx] = 0
            posIdx = posIdx + 1
        #检查循环过程中是否出错
        # print(calRelIdx, len(self.fixRelCpu))
        if calRelIdx != len(self.fixRelCpu):
            print('FixData iteration wrong in updateFixData!')
            print(calRelIdx,len(self.fixRelCpu))
            exit()
        #对所有要进一步修复的断点进行数据准备
        self.stSearchPosCpu = []
        for st in tmpToSearchSt:
            for adjustPos in range(self.checkIt):
                _2Pos = self._1to2(st)
                _2Pos[0] = _2Pos[0] - 1 - self.checkStep * adjustPos
                _2Pos[1] = _2Pos[1] - 1 - self.checkStep * adjustPos
                toSearch1Pos = self._2to1(_2Pos)
                self.stSearchPosCpu.append(toSearch1Pos)
        for ed in tmpToSearchEd:
            for adjustPos in range(self.checkIt):
                _2Pos = self._1to2(ed)
                _2Pos[0] = _2Pos[0] + self.checkLen + self.checkStep * adjustPos
                _2Pos[1] = _2Pos[1] + self.checkLen + self.checkStep * adjustPos
                toSearch1Pos = self._2to1(_2Pos)
                self.stSearchPosCpu.append(toSearch1Pos)
        if len(self.stSearchPosCpu)>0:
            self.stSearchPosCpu = np.array(self.stSearchPosCpu,dtype=np.int32)
            self.stSearchPos = gpuarray.to_gpu(self.stSearchPosCpu)
            self.fixRelCpu = np.zeros(len(self.stSearchPosCpu),dtype=np.int32)
            self.fixRel = gpuarray.to_gpu(self.fixRelCpu)

    def plotSegment(self):
        stringColor = ['#43A047','#D50000', '#FF9100', '#E6DA5A', '#F07CE8', '#FFD67F']
        pltColor = np.array([mcolors.to_rgba(c) for c in stringColor])
        #先画概率
        fig,ax = plt.subplots()
        probdp = self.probRelGpu.get()
        tmpMin = np.min(probdp[probdp > 0])
        probdp[probdp <= 0] = tmpMin / 10
        sns.heatmap(data=np.array(-np.log10(probdp)).reshape(self.blen, self.alen),
                    cbar=False, xticklabels=False, yticklabels=False)
        #再按照从后往前的顺序画每一步迭代结果

        #画Step3
        segData = []
        typeN = []
        for i in range(len(self._3StepSegStPosCpu)):
            stXy = self._1to2(self._3StepSegStPosCpu[i])
            edXy = self._1to2(self._3StepSegEdPosCpu[i])
            segData.append([(stXy[0], stXy[1]), (edXy[0], edXy[1])])
            typeN.append(pltColor[2])
        lc = mc.LineCollection(segData, colors=np.array(typeN),linewidths=8)
        ax.add_collection(lc)

        #画Step2
        # segData = []
        # typeN = []
        # for i in range(len(self._2StepSegStPosCpu)):
        #     stXy = self._1to2(self._2StepSegStPosCpu[i])
        #     edXy = self._1to2(self._2StepSegEdPosCpu[i])
        #     segData.append([(stXy[0],stXy[1]),(edXy[0],edXy[1])])
        #     typeN.append(pltColor[1])
        # lc = mc.LineCollection(segData,colors=np.array(typeN),linewidths=8)
        # ax.add_collection(lc)
        #画Step1
        # segData = []
        # typeN = []
        # for i in range(len(self._1StepSegStPosCpu)):
        #     stXy = self._1to2(self._1StepSegStPosCpu[i])
        #     edXy = self._1to2(self._1StepSegEdPosCpu[i])
        #     segData.append([(stXy[0], stXy[1]), (edXy[0], edXy[1])])
        #     typeN.append(pltColor[0])
        # lc = mc.LineCollection(segData, colors=np.array(typeN),linewidths=8)
        # ax.add_collection(lc)

        fig.set_size_inches(80*(self.alen/max(self.alen,self.blen)), 80*(self.blen/max(self.alen,self.blen)))
        ax.set_ylim([0,self.blen])
        ax.set_xlim([0,self.alen])
        fig.tight_layout(pad=0)
        plt.savefig(f'plot/test200/{self.myId}.png')

        plt.close()

    def mergeOverlapReads(self):

        def solveNowList(nowList,mergeUfs):
            for i in range(len(nowList)):
                j = i+1
                while j<len(nowList):
                    a = nowList[i].st[0]
                    b = nowList[i].ed[0]
                    c = nowList[j].st[0]
                    d = nowList[j].ed[0]
                    if (c>=a and c<=b) or (d>=a and d<=b) or \
                       (a>=c and a<=d) or (b>=c and b<=d):
                        mergeUfs.union(nowList[i].id,nowList[j].id)
                    j = j+1

        mergeUfs = ttUFS(len(self.segStPosCpu))
        minInterceptDiff = 3

        nowSegs = []
        for i in range(len(self.segStPosCpu)):
            nowSegs.append(alignSegment(self._1to2(self.segStPosCpu[i]),self._1to2(self.segEdPosCpu[i]),i))
        nowSegs = sorted(nowSegs,key=attrgetter('intercept'))

        nowIntercept = -1e9
        nowList = []
        for i in range(len(nowSegs)):
            if nowSegs[i].intercept-nowIntercept > minInterceptDiff:
                if len(nowList)>1:
                    solveNowList(nowList,mergeUfs)
                nowList = []
            nowIntercept = nowSegs[i].intercept
            nowList.append(nowSegs[i])
        if len(nowList)>1:
            solveNowList(nowList,mergeUfs)

        for i in range(len(self.segStPosCpu)):
            fa = mergeUfs.find_set(i)

            i2Pos = self._1to2(self.segStPosCpu[i])
            fa2Pos = self._1to2(self.segStPosCpu[fa])
            if fa2Pos[0]>i2Pos[0]:
                self.segStPosCpu[fa] = self.segStPosCpu[i]

            i2Pos = self._1to2(self.segEdPosCpu[i])
            fa2Pos = self._1to2(self.segEdPosCpu[fa])
            if fa2Pos[0]<i2Pos[0]:
                self.segEdPosCpu[fa] = self.segEdPosCpu[i]

        tmpSegStPosCpu = []
        tmpSegEdPosCpu = []
        for i in range(len(self.segStPosCpu)):
            if i==mergeUfs.find_set(i):
                tmpSegEdPosCpu.append(self.segEdPosCpu[i])
                tmpSegStPosCpu.append(self.segStPosCpu[i])
        self.segStPosCpu = np.array(tmpSegStPosCpu,dtype=np.int32)
        self.segEdPosCpu = np.array(tmpSegEdPosCpu,dtype=np.int32)

    def connectSegments(self):
        def overlap(sega,segb):
            a = sega.st[0]
            b = sega.ed[0]
            c = segb.st[0]
            d = segb.ed[0]
            if (c>=a and c<=b) or (d>=a and d<=b) or \
               (a>=c and a<=d) or (b>=c and b<=d):
                return True
            else:
                return False

        def solveNowList(nowList,toTestLinkStId,toTestLinkEdId):
            toTestLinkStId = []
            toTestLinkEdId = []
            for i in range(len(nowList)):
                j = i+1
                while j<len(nowList):
                    if (not overlap(nowList[i],nowList[j])) and \
                        abs(nowList[i].intercept - nowList[j].intercept) <= 50:
                        if nowList[i].st < nowList[j].st:
                            toTestLinkStId.append(nowList[i].id)
                            toTestLinkEdId.append(nowList[j].id)
                        else:
                            toTestLinkStId.append(nowList[j].id)
                            toTestLinkEdId.append(nowList[i].id)
                    j = j+1
            return toTestLinkStId,toTestLinkEdId
        toTestLinkStId = []
        toTestLinkEdId = []
        nowSegs = []
        for i in range(len(self.segStPosCpu)):
            nowSegs.append(alignSegment(self._1to2(self.segStPosCpu[i]), self._1to2(self.segEdPosCpu[i]), i))
        nowSegs = sorted(nowSegs, key=attrgetter('intercept'))

        # for seg in nowSegs:
            # print(seg.id,self._1to2(seg.st)[0],self._1to2(seg.st)[1],self._1to2(seg.ed)[0],self._1to2(seg.ed)[1])
            # print(seg.id, seg.st, seg.ed, seg.intercept)

        minInterceptDiff = 50
        nowIntercept = -1e9

        nowList = []
        for i in range(len(nowSegs)):
            if nowSegs[i].intercept - nowIntercept > minInterceptDiff:
                if len(nowList) > 1:
                    tmpStId,tmpEdId = solveNowList(nowList, toTestLinkStId, toTestLinkEdId)
                    toTestLinkStId = toTestLinkStId + tmpStId
                    toTestLinkEdId = toTestLinkEdId + tmpEdId
                nowList = []
            nowIntercept = nowSegs[i].intercept
            nowList.append(nowSegs[i])
        if len(nowList) > 1:
            tmpStId, tmpEdId = solveNowList(nowList, toTestLinkStId, toTestLinkEdId)
            toTestLinkStId = toTestLinkStId + tmpStId
            toTestLinkEdId = toTestLinkEdId + tmpEdId

        if len(toTestLinkStId)<=1:
            return
        kerCheckPairRela = getCKerFunction(f'{kerP}/dotPlotSvKer.cpp','checkPairRela')
        toTestLinkStId = np.array(toTestLinkStId,dtype=np.int32)
        toTestLinkEdId = np.array(toTestLinkEdId,dtype=np.int32)
        toTestLinkStIdGpu = gpuarray.to_gpu(toTestLinkStId)
        toTestLinkEdIdGpu = gpuarray.to_gpu(toTestLinkEdId)
        self.segStPos = gpuarray.to_gpu(self.segStPosCpu)
        self.segEdPos = gpuarray.to_gpu(self.segEdPosCpu)
        relaRel = np.zeros(len(toTestLinkEdIdGpu),dtype=np.int32)
        relaRelGpu = gpuarray.to_gpu(relaRel)
        N = np.int32(len(relaRel))
        kerCheckPairRela(toTestLinkStIdGpu,toTestLinkEdIdGpu,self.segStPos,self.segEdPos,self.probRelGpu,relaRelGpu,self.alen,self.blen,N,block=(1024,1,1),grid=(int(np.ceil(N/1024)),1))

        myUfs = ttUFS(len(self.segStPosCpu))
        relaRel = relaRelGpu.get()
        for i in range(N):
            if relaRel[i]==1:
                myUfs.union(toTestLinkStId[i],toTestLinkEdId[i])

        for i in range(len(self.segStPosCpu)):
            fa = myUfs.find_set(i)

            i2Pos = self._1to2(self.segStPosCpu[i])
            fa2Pos = self._1to2(self.segStPosCpu[fa])
            if fa2Pos[0]>i2Pos[0]:
                self.segStPosCpu[fa] = self.segStPosCpu[i]

            i2Pos = self._1to2(self.segEdPosCpu[i])
            fa2Pos = self._1to2(self.segEdPosCpu[fa])
            if fa2Pos[0]<i2Pos[0]:
                self.segEdPosCpu[fa] = self.segEdPosCpu[i]

        tmpSegStPosCpu = []
        tmpSegEdPosCpu = []
        for i in range(len(self.segStPosCpu)):
            if i == myUfs.find_set(i):
                tmpSegEdPosCpu.append(self.segEdPosCpu[i])
                tmpSegStPosCpu.append(self.segStPosCpu[i])

        self.segStPosCpu = np.array(tmpSegStPosCpu,dtype=np.int32)
        self.segEdPosCpu = np.array(tmpSegEdPosCpu,dtype=np.int32)

        nowSegs = []
        for i in range(len(self.segStPosCpu)):
            nowSegs.append(alignSegment(self._1to2(self.segStPosCpu[i]), self._1to2(self.segEdPosCpu[i]), i))
        nowSegs = sorted(nowSegs, key=attrgetter('intercept'))
        # for seg in nowSegs:
            # print(seg.id,self._1to2(seg.st)[0],self._1to2(seg.st)[1],self._1to2(seg.ed)[0],self._1to2(seg.ed)[1])
            # print(seg.intercept)
