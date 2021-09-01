//Func generateDp
extern "C"{
__global__ void generateDp(char *seqa, char *seqb, char *rel, int alen, int blen, int N,int width,int minOverlap){
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if(i>=N){
		return;
    }
	int sum = 0;
	int stb = i/alen;
	int sta = i%alen;
	for(int j=stb,k=sta;j<min(int(stb+width),int(blen))&&k<min(int(sta+width),int(alen));j++,k++){
		sum += (seqb[j]==seqa[k]);
	}
	if(sum>=minOverlap){
		rel[i]++;
	}
}
}
//Endfunc

//Func generateProbDp
#include <math.h>
extern "C"{
__global__ void generateProbDp(char *seqa, char *seqb, float *rel, int alen, int blen, int N,int bin){
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if(i>=N){
		return;
    }
	int sum = 0;
	int stb = i/alen;
	int sta = i%alen;
	int validBin = 0;
	float agctNum[4];
	for(int j=0;j<4;j++){
		agctNum[j] = 0;
	}
	for(int j=stb,k=sta;j<min(int(stb+bin),int(blen))&&k<min(int(sta+bin),int(alen));j++,k++){
		sum += (seqb[j]==seqa[k]);
		validBin++;
		switch(seqa[k]){
			case 'A':agctNum[0]++;break;
			case 'G':agctNum[1]++;break;
			case 'C':agctNum[2]++;break;
			case 'T':agctNum[3]++;break;
//			default: printf('Wrong AGCT!\n');
		}
		switch(seqb[j]){
			case 'A':agctNum[0]++;break;
			case 'G':agctNum[1]++;break;
			case 'C':agctNum[2]++;break;
			case 'T':agctNum[3]++;break;
//			default: printf('Wrong AGCT!\n');
		}
	}
	if(validBin<bin){
		rel[i] = 1;
		return;
	}
	for(int j=0;j<4;j++){
		agctNum[j] /= 2*validBin;
	}
	double diffp = 0;
	for(int j=0;j<4;j++){
		diffp += agctNum[j]*(1-agctNum[j]);
	}
	double samep = 1-diffp;
	double sump = 0;
	double cpart = 1;
	double diffPart = pow(diffp,validBin);
	double samePart = 1;
	for(int j=0;j<=sum;j++){
		sump += cpart*diffPart*samePart;
		cpart *= validBin-j;
		cpart /= j+1;
		diffPart /= diffp;
		samePart *= samep;
	}
	rel[i] = 1-sump;
	if(rel[i]<0){
		rel[i] = 0;
	}
}
}
//Endfunc

//Func initStEd
extern "C"{
__global__ void initStEd(float *prob, int *st, int *ed, int N){
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if(i>=N || prob[i]>1e-6){
		return;
    }
    if(prob[i]>1e-6){
		st[i] = -1;
		ed[i] = -1;
    }
    else{
    	st[i] = i;
    	ed[i] = i;
    }
}
}
//Endfunc


//Func findPairRelaForward
extern "C"{
__global__ void findPairRelaForward(int *ttFlag,int *segStPos, int *st, int *ed, int *ned,int alen, int blen, int N){
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if(i>=N){
		return;
	}
	i = int(ed[segStPos[i]]);
	if(ed[i]!=i){
//	    printf("NO mathch!\n");
        return;
	}
	ttFlag[i] = i;
	ned[i] = -1;
	//The prob for length test
	double limitP=0.2;
	int maxLimitSegLen=10;
	int pb = i/alen;
	int pa = i%alen;
	int totLen = alen*blen;
	int stpb = st[i]/alen;
	int stpa = st[i]%alen;
	int myLen = min(abs(pb-stpb),abs(pa-stpa))*1.14;
	myLen += 33;
//	int searchRange = min(15,myLen/10);
	int searchRange = alen/5;
	for(int j=1;j<=searchRange;j++){
		int sideRange = max(1,j/10);
		sideRange = min(3,sideRange);
		for(int k=0;k<=sideRange;k++){
			int spa = pa+j;
			int spb = pb+j+k;
			int si = spb*alen+spa;
			if(si<totLen && st[si]==si){
				double d = j;
				double a = ed[si]/alen - stpb;
				double p = (2/a)*d - (d*d/(a*a));
				if(j<=maxLimitSegLen || p<=limitP){
					ned[i] = si;
					return;
				}
			}
			spb = pb+j-k;
			si = spb*alen+spa;
			if(si<totLen && st[si]==si){
				double d = j;
				double a = ed[si]/alen - stpb;
				double p = (2/a)*d - (d*d/(a*a));
				if(j<=maxLimitSegLen || p<=limitP){
					ned[i] = si;
					return;
				}
			}
		}
	}
}
}
//Endfunc


//Func findPairRelaBackward
extern "C"{
__global__ void findPairRelaBackward(int *ttFlag,int *segStPos, int *st, int *ed, int *nst,int alen, int blen, int N){
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if(i>=N){
		return;
	}
	i=segStPos[i];
	if(st[i]!=i){
//	    printf("No match!!!!!!!!!\n");
        return;
	}
	ttFlag[i] = i;
	nst[i] = -1;
	//The prob for length test
	double limitP=0.2;
	int maxLimitSegLen=10;
	int pb = i/alen;
	int pa = i%alen;
	int edpb = ed[i]/alen;
	int edpa = ed[i]%alen;
	int myLen = min(abs(pb-edpb),abs(pa-edpa))*1.14;
	myLen += 33;
//	int searchRange = min(15,myLen/10);
	int searchRange = alen/5;
	for(int j=1;j<=searchRange;j++){
		int sideRange = max(5,j/10);
		sideRange = min(3,sideRange);
		for(int k=0;k<=sideRange;k++){
			int spa = pa-j;
			int spb = pb-j-k;
			int si = spb*alen+spa;
			if(si>=0 && ed[si]==si){
				double d = j;
				double a = edpb - st[si]/alen;
				double p = (2/a)*d - (d*d/(a*a));
				if(j<=maxLimitSegLen || p<=limitP){
					nst[i] = si;
					return;
				}
			}
			spb = pb-j+k;
			si = spb*alen+spa;
			if(si>=0 && ed[si]==si){
				double d = j;
				double a = edpb - st[si]/alen;
				double p = (2/a)*d - (d*d/(a*a));
				if(j<=maxLimitSegLen || p<=limitP){
					nst[i] = si;
					return;
				}
			}
		}
	}
}
}
//Endfunc

//Func findPairPair
extern "C"{
__global__ void findPairPair(int *nst, int *ned, int *flag, int *st, int *ed,int N){
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if(i>=N){
		return;
	}
	flag[i]=0;
	if(nst[i]>=0 && ned[ nst[i] ] == i){
		st[i] = nst[i];
		flag[i]++;
	}
	if(ned[i]>=0 && nst[ ned[i] ] == i){
		ed[i] = ned[i];
	}
}
}
//Endfunc

//Func updateStEd
extern "C"{
__device__ int lower_bound(int *st,int *ed,int tar){
    int mid, left =0, right = ed - st - 1;
    while(left <= right){
        mid = (left+right) >> 1;
        if(*(st+mid) >= tar) right = mid-1;
        else left = mid+1;
    }
    return left;
}
__global__ void updateStEd(int *st, int *ed, int *segStFlag, int *segStPos,int segNum,int N){
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if(i>=N || st[i]<0){
		return;
	}

	int root = i;
	while(st[root] != root){
		root=st[root];
	}
	st[i] = root;
	int inSegStPos=0;
	int idx = lower_bound(segStPos,segStPos+segNum,i);
    int nowIdx = segNum;

	//if i is in present segStPos
	if(idx<segNum && segStPos[idx]==i){
		//weather i is a starter after the latest iteration or not
		//the segStFlag should be update
		if(st[i]==i){
			segStFlag[idx] = 1;
		}else{
			segStFlag[idx] = 0;
		}
	}
	root = i;
	while(ed[root] != root){
		root=ed[root];
	}
	ed[i] = root;
}
}
//Endfunc

//Func staSegNum
extern "C"{
__global__ void staSegNum(int *st, int *flag, int N){
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if(i>=N){
		return;
	}
	flag[i]	= 0;
	if(st[i]==i){
		flag[i]	= 1;
	}
}
}
//Endfunc

//Func fixCheck
extern "C"{
__global__ void fixCheck(int *pos,float *prob,float *rel,int totLen,int alen,int blen,int N){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i>=N){
        return;
    }
    int np = pos[i];
    if(np<0 || np>=totLen){
        rel[i] = 1;
        return;
    }
    int pa = np%alen;
    int pb = np/alen;
    double inSumP=0,outSumP=0;
    int inNum=0,outNum=0;
    for(int k=-4;k<=0;k++){
        int ctPa = pa+k;
        int ctPb = pb+k;
        if(ctPa<0 || ctPb<0 || ctPa>=alen || ctPb>=blen){
            continue;
        }
        for(int j=-3;j<=3;j++){
            int npa = ctPa+j;
            int npb = ctPb;
            if(npa<0 || npa>=alen){
                continue;
            }
            int si = npa+npb*alen;
            if(j>1 || j<-1){
                outSumP+=prob[si];
                outNum++;
            }else{
                inSumP+=prob[si];
                inNum++;
            }
        }
    }
    if(inNum<3 || outNum<10){
        rel[i] = 1;
        return;
    }


//    rel[i] = (inSumP/inNum)/(outSumP/outNum);
//    return;

    rel[i] = (inSumP/inNum)/(outSumP/outNum);

}
}
//Endfunc

//Func checkPairRela
extern "C"{

__device__ int calInterceptSub(int _1pa,int _1pb,int _2pa,int _2pb){
    int interceptSub = (_1pb-_1pa) - (_2pb-_2pa);
    if(interceptSub<0){
        interceptSub = -interceptSub;
    }
    return interceptSub;
}

__device__ bool overlap(int _1pa,int _2pa){
    return _1pa>=_2pa;
}

__device__ double calDis(int nstPa, int nstPb, int nedPa, int nedPb){
    if(nstPa > nedPa) return 0;
    else return sqrtf((nstPa-nedPa)*(nstPa-nedPa) + (nstPb-nedPb)*(nstPb-nedPb));
}
__global__ void checkPairRela(int *stPointId,int *edPointId,int *segStPos,int *segEdPos,
                                float *prob,int *rel,int alen,int blen,int toCheckNum){
    int ind = blockIdx.x * blockDim.x + threadIdx.x;
    if(ind>=toCheckNum){
        return;
    }
    int stPb = segEdPos[stPointId[ind]]/alen;
    int stPa = segEdPos[stPointId[ind]]%alen;
    int edPb = segStPos[edPointId[ind]]/alen;
    int edPa = segStPos[edPointId[ind]]%alen;
    int nstPb=stPb, nstPa=stPa, nedPb=edPb, nedPa=edPa;
    int initInterceptSub = calInterceptSub(stPa,stPb,edPa,edPb);
    int searchRange = 5;
    int maxSub = max(int(initInterceptSub*1.1),initInterceptSub+searchRange);
    int findNewStPoint,findNewEdPoint;
    while(calInterceptSub(nstPa,nstPb,nedPa,nedPb)<maxSub &&
            (!overlap(nstPa,nedPa)) ){
        //以下计算新的st位点
        double nowMinVal = 1;
        findNewStPoint = 1;
        int tmpPa,tmpPb,tarPa=-1,tarPb=-1,_1Pos;
        for(int i=searchRange;i>0;i--){
            for(int j=1;j<=searchRange;j++){
                tmpPa = nstPa+i;
                tmpPb = nstPb+j;
                if(tmpPa>=0 && tmpPa<alen && tmpPb>=0 && tmpPb<blen){
                    _1Pos = tmpPb*alen+tmpPa;
                    if(prob[_1Pos]<nowMinVal){
                        nowMinVal = prob[_1Pos];
                        tarPa = tmpPa;
                        tarPb = tmpPb;
                    }
                }
            }
        }
        if(tarPa!=-1 && nowMinVal<1e-2){
            nstPa = tarPa;
            nstPb = tarPb;
        }else{
            findNewStPoint = 0;
        }

        //以下计算新的ed位点
        nowMinVal = 1;
        findNewEdPoint = 1;
        tarPa=-1;
        tarPb=-1;
        for(int i=searchRange;i>0;i--){
            for(int j=1;j<=searchRange;j++){
                tmpPa = nedPa-i;
                tmpPb = nedPb-j;
                if(tmpPa>=0 && tmpPa<alen && tmpPb>=0 && tmpPb<blen){
                    _1Pos = tmpPb*alen+tmpPa;
                    if(prob[_1Pos]<nowMinVal){
                        nowMinVal = prob[_1Pos];
                        tarPa = tmpPa;
                        tarPb = tmpPb;
                    }
                }
            }
        }
        if(tarPa!=-1 && nowMinVal<1e-2){
            nedPa = tarPa;
            nedPb = tarPb;
        }else{
            findNewEdPoint = 0;
        }

        if((!findNewStPoint) && (!findNewEdPoint)){
            break;
        }
    }
//    if(stPointId[ind]==123 || edPointId[ind]==123){
//        printf("==================================\n");
//        printf("%d\t%d\t%d\t%d\n",nstPa,nstPb,nedPa,nedPb);
//        printf("%d\n",calInterceptSub(nstPa,nstPb,nedPa,nedPb));
//        printf("%.4lf\n",calDis(nstPa,nstPb,nedPa,nedPb));
//        printf("==================================\n");
//    }
    if(calInterceptSub(nstPa,nstPb,nedPa,nedPb)<=searchRange &&
       calDis(nstPa,nstPb,nedPa,nedPb)<=searchRange){
//        printf("%d\t%d\t%d\t%d\n",nstPa,nstPb,nedPa,nedPb);
//        printf("%.4lf\t%.4lf\n",calInterceptSub(nstPa,nstPb,nedPa,nedPb),calDis(nstPa,nstPb,nedPa,nedPb));
        rel[ind] = 1;
    }else{
        rel[ind] = 0;
    }

}
}
//Endfunc
