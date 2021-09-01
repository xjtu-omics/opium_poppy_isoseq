from Bio import Entrez
import os
import re
import sys
from ttlib.netBase import myFtp
Entrez.email = "xu_tun@qq.com"

def getLineageInfo(orgnism):
    handle = Entrez.esearch(db="taxonomy", term=orgnism)
    rel = Entrez.read(handle)
    taxId = rel['IdList'][0]
    handle = Entrez.efetch(db='taxonomy',id=taxId,retmode='xml')
    rel = Entrez.read(handle)
    return(rel[0]['Lineage'])

def testNcbiHaveIsoReads(organism):
    searchTerm = '((' + organism + r') AND "transcriptomic"[Source]) AND "pacbio smrt"[Platform]'
    handle = Entrez.esearch(db="sra", term=searchTerm)
    record = Entrez.read(handle)
    p = re.compile(r"'RetMax':\s'[0-9]+'")
    tmp = p.findall(str(record))
    return(tmp[0][11:-1])

def testRefAssemblyIdDict(organismList,outList):
    with open(outList,'w'):
        for org in organismList:
            handle = Entrez.esearch(db='assembly',term=org)
            rel = Entrez.read(handle)
            print(rel['IdList'])

#The assemblyId should be the grey tab named "IDs:\s([0-9]+?)\s"
def downloadReference(assemblyId,dwDir):
    handle = Entrez.esummary(db='assembly',id=assemblyId)
    rel = Entrez.read(handle)
    dirPath = re.findall("'FtpPath_GenBank':\s'.*?nih\.gov/(.*?)'",str(rel))[0]
    ncbiHost = 'ftp.ncbi.nlm.nih.gov'
    ftp = myFtp(ncbiHost,'anonymous','password')
    ftp.downloadFileTree(dwDir,dirPath)
    ftp.close()

def getAssemblyIdFromName(assemblyName):
    handle = Entrez.esearch(db='assembly', term=assemblyName)
    rel = Entrez.read(handle)
    return(rel['IdList'])

def downloadSra(sraId,dwDir):
    sp = '/data/home/xutun/miniconda3/envs/tt/bin/'
    cmd = sp + 'prefetch -O ' + dwDir + ' --max-size 200GB ' + sraId
    os.system(cmd)