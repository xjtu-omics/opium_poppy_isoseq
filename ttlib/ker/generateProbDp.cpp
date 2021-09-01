#include <math.h>
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