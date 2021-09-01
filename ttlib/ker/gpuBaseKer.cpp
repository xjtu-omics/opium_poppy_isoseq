//Func getSumKer
extern "C"{
__global__ void getSumKer(int *tmp,int it,int N){
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if(i>=N){
		return;
	}
	int mergeIdx = (i|(1<<(it-1)));
	int id=1;
	while(mergeIdx>=N && id<=it){
		id++;
		mergeIdx = (i|(1<<(it-id)));
	}
	if((i&((1<<it)-1)) == 0){
		if(mergeIdx < N){
			tmp[i] += tmp[mergeIdx];
		}
	}
}
}
//Endfunc

//Func initArray
extern "C"{
__global__ void initArray(int *tmp,int initValue,int N){
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if(i>=N){
		return;
	}
	tmp[i] = initValue;
}
}
//Endfunc
