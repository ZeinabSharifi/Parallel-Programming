//ONLY MODIFY THIS FILE!
//YOU CAN MODIFY EVERYTHING IN THIS FILE!

#include "fft.h"

#define tx threadIdx.x
#define ty threadIdx.y
#define tz threadIdx.z

#define bx blockIdx.x
#define by blockIdx.y
#define bz blockIdx.z

#define TILE 1024
//#define TILE 4

// you may define other parameters here!
// you may define other macros here!
// you may define other functions here!

__global__ void kernelReverseBit2(float* x_r_d, float* x_i_d, const unsigned int M){
	int i = (TILE/2)*bx + tx;
	float tmp_r,tmp_i;
	unsigned int count ;
	unsigned int index, tmp_index ;
	unsigned int reverse_index;
	for(int j=0; j<8; j++){
		index = 8*i+j;
		tmp_index = index;
		count = M;
		reverse_index = 0;
		while(tmp_index)
		{
			reverse_index <<= 1; 
			reverse_index |= tmp_index & 1;
			tmp_index >>= 1;
			count--;
		}
		reverse_index <<= count;

		//printf("index is: %d and reverse_index is %d\n", 8*i+j, reverse_index);

		if(reverse_index > index){
			tmp_r = x_r_d[index];
			tmp_i = x_i_d[index];
			x_r_d[index] = x_r_d[reverse_index];
			x_i_d[index] = x_i_d[reverse_index];
			x_r_d[reverse_index] = tmp_r;
			x_i_d[reverse_index] = tmp_i;
		}
	}

}

__global__ void kernelReverseBit4(float* x_r_d, float* x_i_d, const unsigned int M){
	int i = (TILE/2)*bx + tx;
	float tmp_r,tmp_i;
	unsigned int count ;
	unsigned int index, tmp_index ;
	unsigned int reverse_index;
	for(int j=0; j<8; j++){
		index = 8*i+j;
		tmp_index = index;
		count = M/2;
		reverse_index = 0;
		while(tmp_index)
		{
			reverse_index <<= 2; 
			reverse_index |= tmp_index & 3;
			tmp_index >>= 2;
			count--;
		}
		reverse_index <<= 2*count;

		//printf("index is: %d and reverse_index is %d\n", 8*i+j, reverse_index);

		if(reverse_index > index){
			tmp_r = x_r_d[index];
			tmp_i = x_i_d[index];
			x_r_d[index] = x_r_d[reverse_index];
			x_i_d[index] = x_i_d[reverse_index];
			x_r_d[reverse_index] = tmp_r;
			x_i_d[reverse_index] = tmp_i;
		}
	}
}

__global__ void kernelFFT2Shared(float* x_r_d, float* x_i_d){
	__shared__ float x_r_s[TILE];
	__shared__ float x_i_s[TILE];

	int i = TILE*bx + tx;

	int index;
	float angle;
	float w_r,w_i,tmp_r1,tmp_r2,tmp_i1,tmp_i2;

	x_r_s[tx] = x_r_d[i];
	x_r_s[tx+TILE/2] = x_r_d[i+TILE/2];
	x_i_s[tx] = x_i_d[i];
	x_i_s[tx+TILE/2] = x_i_d[i+TILE/2];



	for(int j=1; j<TILE; j*=2){
		__syncthreads();
		/*
		if(j==1){
			printf("x_r_s[%d] %f x_r_s[%d] %f\n",tx,x_r_s[tx],tx+TILE/2,x_r_s[tx+TILE/2]);
		}
		
		*/
		index = 2*j*(tx/j) + (tx%j) ;
		angle = index%(2*j)*((float)1/(2*j));
		w_r = cos(2*PI*angle);
		w_i = -sin(2*PI*angle);
		//printf("w_r %f w_i %f\n",w_r,w_i);

		//printf("step %d thread %d index %d angle %f\n",j,i,index,angle );

		tmp_r1 = x_r_s[index];
		tmp_i1 = x_i_s[index];
		tmp_r2 = x_r_s[index+j]*w_r -  x_i_s[index+j]*w_i;
		tmp_i2 = x_r_s[index+j]*w_i +  x_i_s[index+j]*w_r;
		//printf("tmp_r1 %f tmp_r2 %f\n",tmp_r1,tmp_r2);

		x_r_s[index] = tmp_r1 + tmp_r2;
		x_i_s[index] = tmp_i1 + tmp_i2;
		x_r_s[index+j] = tmp_r1 - tmp_r2;
		x_i_s[index+j] = tmp_i1 - tmp_i2;
		//printf("x_r_s[%d] %f x_r_s[%d] %f\n",index,x_r_s[index],index+j,x_r_s[index+j]);

	}
	__syncthreads();

	x_r_d[i] = x_r_s[tx];
	x_i_d[i] = x_i_s[tx];
	x_r_d[i+TILE/2] = x_r_s[tx+TILE/2];
	x_i_d[i+TILE/2] = x_i_s[tx+TILE/2];
}

__global__ void kernelFFT2Global(float* x_r_d, float* x_i_d, const unsigned int j){


	int i = bx*(TILE/2)+ tx;

	int index;
	float angle;
	float w_r, w_i,tmp_r1,tmp_r2,tmp_i1,tmp_i2;


	index = (i%j) + (i/j)*j*2 ;
	//printf("thread %d index %d\n",i,index);
	angle = index%(2*j)*((float)1/(2*j));
	w_r = cos(2*PI*angle);
	w_i = -sin(2*PI*angle);
	//printf("thread %d index %d\n");

	tmp_r1 = x_r_d[index];
	tmp_i1 = x_i_d[index];
	tmp_r2 = x_r_d[index+j]*w_r -  x_i_d[index+j]*w_i;
	tmp_i2 = x_r_d[index+j]*w_i +  x_i_d[index+j]*w_r;

	x_r_d[index] = tmp_r1 + tmp_r2;
	x_i_d[index] = tmp_i1 + tmp_i2;
	x_r_d[index+j] = tmp_r1 - tmp_r2;
	x_i_d[index+j] = tmp_i1 - tmp_i2;

}

__global__ void kernelFFT4Shared(float* x_r_d, float* x_i_d){
	__shared__ float x_r_s[TILE];
	__shared__ float x_i_s[TILE];

	int i = by*gridDim.x*TILE + bx*TILE + tx;
	//int i = bx*TILE + tx;


	int index;
	float angle;
	float w_r1,w_i1, w_r2,w_i2, w_r3,w_i3;
	float tmp_r1,tmp_i1,tmp_r2,tmp_i2, tmp_r3,tmp_i3,tmp_r4,tmp_i4;

	x_r_s[tx] = x_r_d[i];
	x_i_s[tx] = x_i_d[i];
	x_r_s[tx+TILE/4] = x_r_d[i+TILE/4];
	x_i_s[tx+TILE/4] = x_i_d[i+TILE/4];
	x_r_s[tx+TILE/2] = x_r_d[i+TILE/2];
	x_i_s[tx+TILE/2] = x_i_d[i+TILE/2];
	x_r_s[tx+3*TILE/4] = x_r_d[i+3*TILE/4];
	x_i_s[tx+3*TILE/4] = x_i_d[i+3*TILE/4];



	for(int j=1; j<TILE; j*=4){
		__syncthreads();
		/*
		if(j==1){
			printf("x_r_s[%d] %f x_r_s[%d] %f\n",tx,x_r_s[tx],tx+TILE/2,x_r_s[tx+TILE/2]);
		}
		
		*/
		index = 4*j*(tx/j) + (tx%j) ;
		angle = index%(4*j)*((float)1/(4*j));

		w_r1 = cos(2*PI*angle);
		w_i1 = -sin(2*PI*angle);
		w_r2 = cos(4*PI*angle);
		w_i2 = -sin(4*PI*angle);
		w_r3 = cos(6*PI*angle);
		w_i3 = -sin(6*PI*angle);

		//printf("w_r %f w_i %f\n",w_r,w_i);

		//printf("step %d thread %d index %d angle %f\n",j,i,index,angle );

		tmp_r1 = x_r_s[index];
		tmp_i1 = x_i_s[index];
		tmp_r2 = x_r_s[index+j]*w_r1 -  x_i_s[index+j]*w_i1;
		tmp_i2 = x_r_s[index+j]*w_i1 +  x_i_s[index+j]*w_r1;
		tmp_r3 = x_r_s[index+2*j]*w_r2 -  x_i_s[index+2*j]*w_i2;
		tmp_i3 = x_r_s[index+2*j]*w_i2 +  x_i_s[index+2*j]*w_r2;
		tmp_r4 = x_r_s[index+3*j]*w_r3 -  x_i_s[index+3*j]*w_i3;
		tmp_i4 = x_r_s[index+3*j]*w_i3 +  x_i_s[index+3*j]*w_r3;
		//printf("tmp_r1 %f tmp_r2 %f\n",tmp_r1,tmp_r2);

		x_r_s[index] = tmp_r1 + tmp_r2 + tmp_r3 + tmp_r4;
		x_i_s[index] = tmp_i1 + tmp_i2 + tmp_i3 + tmp_i4;

		x_r_s[index+j] = tmp_r1 + tmp_i2 - tmp_r3 - tmp_i4;
		x_i_s[index+j] = tmp_i1 - tmp_r2 - tmp_i3 + tmp_r4;

		x_r_s[index+2*j] = tmp_r1 - tmp_r2 + tmp_r3 - tmp_r4;
		x_i_s[index+2*j] = tmp_i1 - tmp_i2 + tmp_i3 - tmp_i4;

		x_r_s[index+3*j] = tmp_r1 - tmp_i2 - tmp_r3 + tmp_i4;
		x_i_s[index+3*j] = tmp_i1 + tmp_r2 - tmp_i3 - tmp_r4;		
		//printf("x_r_s[%d] %f x_r_s[%d] %f\n",index,x_r_s[index],index+j,x_r_s[index+j]);

	}
	__syncthreads();

	x_r_d[i] = x_r_s[tx];
	x_i_d[i] = x_i_s[tx];
	x_r_d[i+TILE/4] = x_r_s[tx+TILE/4];
	x_i_d[i+TILE/4] = x_i_s[tx+TILE/4];
	x_r_d[i+TILE/2] = x_r_s[tx+TILE/2];
	x_i_d[i+TILE/2] = x_i_s[tx+TILE/2];	
	x_r_d[i+3*TILE/4] = x_r_s[tx+3*TILE/4];
	x_i_d[i+3*TILE/4] = x_i_s[tx+3*TILE/4];
}

__global__ void kernelFFT4Global(float* x_r_d, float* x_i_d, const unsigned int j){


	int i = bx*(TILE/2) + tx;

	int index;
	float angle;
	float w_r1,w_i1, w_r2,w_i2, w_r3,w_i3;
	float tmp_r1,tmp_i1,tmp_r2,tmp_i2, tmp_r3,tmp_i3,tmp_r4,tmp_i4;


	index = 4*j*(i/j) + (i%j) ;
	angle = index%(4*j)*((float)1/(4*j));
	//printf("step %d thread %d index %d angle %f\n",j,i,index,angle );

	w_r1 = cos(2*PI*angle);
	w_i1 = -sin(2*PI*angle);
	w_r2 = cos(4*PI*angle);
	w_i2 = -sin(4*PI*angle);
	w_r3 = cos(6*PI*angle);
	w_i3 = -sin(6*PI*angle);
	//printf("i %d w_r1 %f w_i1 %f\n",i,w_r1,w_i1);
	//printf("i %d w_r2 %f w_i2 %f\n",i,w_r2,w_i2);
	//printf("i %d w_r3 %f w_i3 %f\n",i,w_r3,w_i3);

	

	tmp_r1 = x_r_d[index];
	tmp_i1 = x_i_d[index];
	tmp_r2 = x_r_d[index+j]*w_r1 -  x_i_d[index+j]*w_i1;
	tmp_i2 = x_r_d[index+j]*w_i1 +  x_i_d[index+j]*w_r1;
	tmp_r3 = x_r_d[index+2*j]*w_r2 -  x_i_d[index+2*j]*w_i2;
	tmp_i3 = x_r_d[index+2*j]*w_i2 +  x_i_d[index+2*j]*w_r2;
	tmp_r4 = x_r_d[index+3*j]*w_r3 -  x_i_d[index+3*j]*w_i3;
	tmp_i4 = x_r_d[index+3*j]*w_i3 +  x_i_d[index+3*j]*w_r3;

	//printf("i %d tmp_r1 %f tmp_i1 %f\n",i,tmp_r1,tmp_i1);
	//printf("i %d tmp_r2 %f tmp_i2 %f\n",i,tmp_r2,tmp_i2);
	//printf("i %d tmp_r3 %f tmp_i3 %f\n",i,tmp_r3,tmp_i3);
	//printf("i %d tmp_r4 %f tmp_i4 %f\n",i,tmp_r4,tmp_i4);


	x_r_d[index] = tmp_r1 + tmp_r2 + tmp_r3 + tmp_r4;
	x_i_d[index] = tmp_i1 + tmp_i2 + tmp_i3 + tmp_i4;

	x_r_d[index+j] = tmp_r1 + tmp_i2 - tmp_r3 - tmp_i4;
	x_i_d[index+j] = tmp_i1 - tmp_r2 - tmp_i3 + tmp_r4;

	x_r_d[index+2*j] = tmp_r1 - tmp_r2 + tmp_r3 - tmp_r4;
	x_i_d[index+2*j] = tmp_i1 - tmp_i2 + tmp_i3 - tmp_i4;

	x_r_d[index+3*j] = tmp_r1 - tmp_i2 - tmp_r3 + tmp_i4;
	x_i_d[index+3*j] = tmp_i1 + tmp_r2 - tmp_i3 - tmp_r4;		

	//printf("x_r_d[%d] %f x_i_d[%d] %f\n",index,x_r_d[index],index,x_i_d[index]);
	//printf("x_r_d[%d] %f x_i_d[%d] %f\n",index+j,x_r_d[index+j],index+j,x_i_d[index+j]);
	//printf("x_r_d[%d] %f x_i_d[%d] %f\n",index+2*j,x_r_d[index+2*j],index+2*j,x_i_d[index+2*j]);
	//printf("x_r_d[%d] %f x_i_d[%d] %f\n",index+3*j,x_r_d[index+3*j],index+3*j,x_i_d[index+3*j]);	


}



//-----------------------------------------------------------------------------
void gpuKernel(float* x_r_d, float* x_i_d, const unsigned int N, const unsigned int M)
{
	// In this function, both inputs and outputs are on GPU.
	// No need for cudaMalloc, cudaMemcpy or cudaFree.

	if(M%2 == 1){ // M = 23 or 25
		kernelReverseBit2<<< N/(4*TILE), TILE/2 >>>(x_r_d, x_i_d, M); // each thread reverses 8 indices
		
		kernelFFT2Shared<<<N/TILE,TILE/2>>>(x_r_d, x_i_d);
		
		for(int j = TILE; j < N; j*=2){
			kernelFFT2Global<<<N/TILE,TILE/2>>>(x_r_d, x_i_d, j);
		}
	}else{ // M = 24 or 26

		kernelReverseBit4<<< N/(4*TILE), TILE/2 >>>(x_r_d, x_i_d, M);// each thread reverses 8 indices

		dim3 dimGrid(N/(2*TILE),2);
		kernelFFT4Shared<<<dimGrid,TILE/4>>>(x_r_d, x_i_d);

		for(int j = TILE; j < N; j*=4){
			kernelFFT4Global<<<N/(2*TILE),TILE/2>>>(x_r_d, x_i_d, j);
		}
		
	}
	

}
