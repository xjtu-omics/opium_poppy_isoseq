import re
def reforMatCollinearity(colF,outF):
    with open(outF,'w') as of:
        print('GA','GB','group',file=of,sep='\t')
        nowGroup = ''
        for line in open(colF,'r'):
            line = line.strip()
            if line[0] == '#':
                if line[0:12] == '## Alignment':
                    nowGroup=re.findall(r"Alignment\s(.+?):",line)[0]
            else:
                line = line.split('\t')
                print(line[1],line[2],nowGroup,file=of,sep='\t')

