import os

class ttWriter:

    def __init__(self,fileName,cacheNum=100000,initTyp = 'w'):
        self.fn = fileName
        if initTyp == 'w' and os.path.exists(self.fn):
            os.remove(self.fn)
        self.cache = []
        self.cacheNum = cacheNum
        self.inCacheNum = 0

    def writ(self,toWriteStr):
        self.cache.append(toWriteStr)
        self.inCacheNum = self.inCacheNum + 1
        if self.inCacheNum >= self.cacheNum:
            with open(self.fn,'a') as of:
                for outStr in self.cache:
                    print(outStr,file=of)
            self.cache = []
            self.inCacheNum = 0

    def close(self):
        if self.inCacheNum > 0:
            with open(self.fn,'a') as of:
                for outStr in self.cache:
                    print(outStr,file=of)
            self.cache = []
            self.inCacheNum = 0
