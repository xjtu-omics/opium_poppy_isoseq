import os
from collections import defaultdict
with open('../data/total_RI_strict.ioe','r') as f,\
     open('../data/ri.gff','w') as of:
    line=f.readline()
    while True:
        line=f.readline()
        if line=='':
            break
        line=line.strip()
        cpline=line
        line=line.split('\t')
        chrm=line[0]
        int_chr=int(chrm.split('r')[1])
        if int_chr>=12:
            chrm=str('unplaced-scaffold_' + str(int_chr-12))
        gene=line[1]
        info=line[2].split(':')[3].split('-')
        strand=line[2].split(':')[5]
        st=int(info[0])+1
        ed=int(info[1])-1
        print(chrm,'dexseq_prepare_annotation.py\texonic_part',st,ed,'.',strand,'.',gene+':'+chrm+':'+str(st)+':'+str(ed),sep='\t',file=of)
os.system("sort -u -k9,9 ../data/ri.gff | sort -k1,1 -k4,4n  >../data/ri.sorted.gff && rm ../data/ri.gff")
tissues=['S1','S2','S3','S4','S5','S6','S7']
for tis in tissues:
    os.system("bamToBed -i ../../2nd_read/lqs/sdata/%s.bam -cigar | /home/xutun/software/filterSamFile/bin/bedToJunction - | awk 'BEGIN{OFS=\"\\t\"}{print $1,$2,$3,$4,$5,$6}' | sort -k1,1 -k2,2n >../data/%s.intron.bed "%(tis,tis))
    os.system("intersectBed -wao -sorted -f 1.0 -a ../data/ri.sorted.gff -b ../data/%s.intron.bed | awk -v OFS=\"\\t\" '{$16<=0?s[$9]+=0:s[$9]+=$14;}END{for (i in s){print i,s[i]}}' | sort -k1 >../data/%s.exonic_parts.exclusion"%(tis,tis)) 
    os.system("bedtools bamtobed -split -i ../../2nd_read/lqs/sdata/%s.bam |sort -k1,1 -k2,2n >../data/%s.bed"%(tis,tis))
    os.system("awk  -v OFS=\"\\t\" '{print $1,$4,$5,$9;}' ../data/ri.sorted.gff |coverageBed -sorted -split -b ../data/%s.bed -a stdin | awk 'BEGIN{OFS=\"\\t\"}{print $1,$2,$3,$3-$2+1,$5,$4}' | sort -k 6 >../data/%s.exonic_parts.inclusion"%(tis,tis))
    os.system("paste ../data/%s.exonic_parts.inclusion ../data/%s.exonic_parts.exclusion | awk -v \"len=150\" 'BEGIN{OFS=\"\\t\";print \"exon_ID\",\"length\",\"inclusion\",\"exclusion\",\"PSI\"}{NIR=$5/($4+len-1);NER=$8/(len-1);print $6,$4,$5,$8,(NIR+NER<=0)?\"NA\":NIR/(NIR+NER);}'>../data/%s.exonic_parts.psi"%(tis,tis,tis))
