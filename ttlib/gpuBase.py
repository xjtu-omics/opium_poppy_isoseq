from pycuda.compiler import SourceModule
import re
import numpy as np
ttlibPath = '/data/home/xutun/mySrc/ttlib'

def getCKerFunction(cf,funcName):
    funcContent = ''
    inFunc = False
    for line in open(cf,'r'):
        if not inFunc:
            matchList = re.findall('^//Func\s(.+)?$',line)
            if len(matchList)>0 and matchList[0] == funcName:
                inFunc = True
        else:
            matchList = re.findall('^//Endfunc$',line)
            if len(matchList) > 0:
                inFunc=False
        if inFunc:
            funcContent = funcContent + line
    funcRel = SourceModule(funcContent,no_extern_c=True,include_dirs=['/data/home/xutun/miniconda3/envs/gpu/include/boost/numeric/odeint/external/thrust/']).get_function(funcName)
    return funcRel

def sumGpu(intArrGpu,N):
    nowIt = 1
    nowRange = 2
    sumKer = getCKerFunction(f'{ttlibPath}/ker/gpuBaseKer.cpp','getSumKer')
    while True:
        sumKer(intArrGpu,np.int32(nowIt),N,block=(1024,1,1),grid=(int(np.ceil(N/1024)),1))
        if nowRange >= N:
            break
        nowIt = nowIt + 1
        nowRange = nowRange * 2
    return intArrGpu.get()[0]
