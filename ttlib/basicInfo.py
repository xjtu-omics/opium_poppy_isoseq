class bainfo:
    def __init__(self,_chr,_st=-1,_ed=-1,_trans='',_strand='+'):
        if ':' not in _chr:
            self.chr=_chr
            self.st=_st
            self.ed=_ed
            self.trans=_trans
            self.strand=_strand
        else:
            self.chr = _chr.split(':')[0]
            self.st = int(_chr.split(':')[1].split('-')[0])
            self.ed = int(_chr.split(':')[1].split('-')[1].split('(')[0])
            self.strand = _chr.split('(')[1][0]
            self.trans = 'nope'

    def __eq__(self,other):
        return self.chr==other.chr and self.st==other.st

    def __lt__(self,other):
        if self.chr!=other.chr:
            return self.chr<other.chr
        else:
            return self.st<other.st
    
    def __cmp__(self,other):
        if self.chr<other.chr:
            return -1
        elif self.chr>other.chr:
            return 1
        elif self.st<other.st:
            return -1
        elif self.st>other.st:
            return 1
        elif self.st==other.st:
            return 0
        else:
            return

    def print(self,of=None):
        if of == None:
            print(self.chr,self.st,self.ed,self.trans,'.',self.strand,sep='\t')
        else:
            print(self.chr,self.st,self.ed,self.trans,'.',self.strand,sep='\t',file=of)

    def getStrBainfo(self):
        return f'{self.chr}:{self.st}-{self.ed}({self.strand})'

class ttSeq:
    def __init__(self,_qn,_seq):
        self.qn = _qn
        self.seq = _seq

def overlapMax(a,b):
    if a.chr != b.chr:
        return 0
    overDis=max(min(a.ed,b.ed)-max(a.st,b.st),0)
    return max(overDis/(a.ed-a.st+1),overDis/(b.ed-b.st+1))

def overlapMin(a,b):
	overDis=max(min(a.ed,b.ed)-max(a.st,b.st),0)
	return min(overDis/(a.ed-a.st+1),overDis/(b.ed-b.st+1))

def overlapA(a,b):
	overDis=max(min(a.ed,b.ed)-max(a.st,b.st),0)
	return overDis/(a.ed-a.st+1)

def mergeBainfo(a,b,checkStrand=False):
    if a.chr != b.chr or overlapMax(a,b) == 0 or (checkStrand and a.strand != b.strand):
        print(f'Merge false of {a.getStrBainfo()} and {b.getStrBainfo()}')
        exit()
    else:
        nst = min(a.st,b.st)
        ned = max(a.ed,b.ed)
        return bainfo(a.chr,nst,ned,a.trans,a.strand)

def calDis(a,b):
    if a.chr != b.chr:
        print(a.trans,b.trans,'not same chr!')
        exit()
    c = [a,b]
    c.sort()
    return c[1].st-c[0].ed
