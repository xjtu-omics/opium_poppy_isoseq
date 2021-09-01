from stringBase import printTabList
def gtfName2Chr(inGtfFileName,outGtfFileName,originMarkPre,baseChrNum):
    with open(outGtfFileName,'w') as of:
        print('# gtfName2Chr(%s,%s,%s,%s)'%(inGtfFileName,outGtfFileName,originMarkPre,baseChrNum),file=of)
        for line in open(inGtfFileName,'r'):
            line = line.strip()
            if line[0] == '#':
                print(line,file=of)
                continue
            ll = line.split('\t')
            if originMarkPre in ll[0]:
                ll[0] = 'chr' + str(int(ll[0].split(originMarkPre)[1]) + baseChrNum + 1)
            printTabList(ll,of)
