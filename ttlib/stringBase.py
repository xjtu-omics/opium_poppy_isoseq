import uuid

def generateWorkD():
    return uuid.uuid4().hex
def generateRandomString():
    return uuid.uuid4().hex

def getReverseComplete(mySeq):
    pos = len(mySeq) - 1
    retStr = ''
    while pos >=0:
        if mySeq[pos] == 'A':
            retStr = retStr + 'T'
        elif mySeq[pos] == 'a':
            retStr = retStr + 't'
        elif mySeq[pos] == 'T':
            retStr = retStr + 'A'
        elif mySeq[pos] == 't':
            retStr = retStr + 'a'
        elif mySeq[pos] == 'C':
            retStr = retStr + 'G'
        elif mySeq[pos] == 'c':
            retStr = retStr + 'g'
        elif mySeq[pos] == 'G':
            retStr = retStr + 'C'
        elif mySeq[pos] == 'g':
            retStr = retStr + 'c'
        else:
            retStr = retStr + 'N'
        pos = pos - 1
    return retStr

def printTabList(myList,of=None):
    ll = len(myList)
    if of == None:
        for i in range(ll):
            print(myList[i],end='')
            if i != ll-1:
                print('\t',end='')
            else:
                print()
    else:
        for i in range(ll):
            print(myList[i],end='',file=of)
            if i != ll-1:
                print('\t',end='',file=of)
            else:
                print('',file=of)

def removePolyAT(seq):
    window = 10
    nst = -1
    cutP = 0.8
    ned = -1
    seqLen = len(seq)
    atNum = 0
    for i in range(seqLen):
        if seq[i]=='A' or seq[i]=='T' or \
           seq[i]=='a' or seq[i]=='t':
            atNum = atNum+1
        if i>=window:
            np = i-window
            if seq[np]=='A' or seq[np]=='T' or \
                seq[np]=='a' or seq[np]=='t':
                atNum = atNum-1
            if atNum < window*cutP:
                nst = i
                break
    if nst==-1:
        return None

    i = seqLen-1
    j = 0
    atNum = 0
    while i>=0:
        if seq[i]=='A' or seq[i]=='T' or \
           seq[i]=='a' or seq[i]=='t':
            atNum = atNum+1
        if j>=window:
            np = i+window
            if seq[np] == 'A' or seq[np] == 'T' or \
               seq[np] == 'a' or seq[np] == 't':
                atNum = atNum - 1
            if atNum < window * cutP:
                ned = i+1
                break
        j = j+1
        i = i-1
    if ned==-1:
        return None

    return seq[nst:ned]

def getGc(seq):
    gcNum = 0
    seqLen = len(seq)
    for i in range(seqLen):
        if seq[i]=='G' or seq[i]=='C' or \
           seq[i]=='g' or seq[i]=='c':
            gcNum = gcNum+1
    return gcNum/seqLen
