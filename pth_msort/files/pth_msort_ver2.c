// Include your C header files here
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "pth_msort.h"

#define CUTOFF 8

int* a;
int* tmp;
unsigned int n;
unsigned int samples[4] = {0}; // 2 samples from each half
unsigned int boudary1[5] = {0}; //boundaries of independant merges in first half
unsigned int boudary2[5] = {0}; //boundaries of independant merges in second half

struct Bound{
	unsigned int lo1;
	unsigned int hi1;
	unsigned int lo2;
	unsigned int hi2;
	unsigned int offset;
	unsigned int length;
};

struct Search{
	unsigned int lo;
	unsigned int hi;
	unsigned int rank;
	int num;
};


void insertionSort(int* thread_tmp, unsigned int lo, unsigned int hi) {
	//printf("lo & hi in insertionSort is : lo %u hi %u\n", lo, hi);
	//printf("lo: %d\n", thread_tmp[lo]);
	//printf("hi: %d\n", thread_tmp[hi]);

	unsigned int k;
	int j;
	unsigned int i;
	for (k = lo+1; k <= hi; k++) {
		int tmp = thread_tmp[k];
		//printf("tmp is%d\n", tmp);
		for (j = k-1; j >= lo && tmp < thread_tmp[j]; j--) {
			printf("thread_tmp[j] is%d\n", thread_tmp[j]);
			printf("thread_tmp[j+1] is%d\n", thread_tmp[j+1]);
			thread_tmp[j+1] = thread_tmp[j];
		}
		thread_tmp[j+1] = tmp;
		//for(i=lo;i<=hi;i++)
		//	printf("thread_tmp %u th element in %u iteration is: %d", i,k-1,thread_tmp[i]);
	}

}

void merge(int* thread_a, int* thread_tmp, unsigned int lo, unsigned int mid, unsigned int hi){
	unsigned int i, k, j;
	/*
	for(k=lo; k<=hi; k++){
		tmp[k] = a[k];
	}
	*/
	i = lo;
	j = mid+1;
	for (k = lo; k <= hi; k++) {
		if (i > mid) {
			thread_tmp[k] = thread_a[j];
			j++;
		}
		else if (j > hi) {
			thread_tmp[k] = thread_a[i];
			i++;
		}
		else if (thread_a[j] < thread_a[i]) {
			thread_tmp[k] = thread_a[j];
			j++;
		}
		else {
			thread_tmp[k] = thread_a[i];
			i++;
		}
	}
}

void mergeSort(int* thread_a, int* thread_tmp, unsigned int lo, unsigned int hi) {

	/* insertion sort did not work
	if (hi - lo <= CUTOFF - 1) {// for subarrays less than 8 use insertion sort instead of recursive calls
		insertionSort(thread_tmp, lo, hi); // like merge sort applied, insertion sort also sorts in tmp array
		unsigned int k;
		for (k = lo; k <= hi; k++)
			printf("local tmp after sort is: %d\n", thread_tmp[k]);
		return;
	}
	*/
	if(hi <= lo) return;
	unsigned int mid = lo + (hi - lo) / 2;

	mergeSort(thread_tmp, thread_a, lo, mid);
	mergeSort(thread_tmp, thread_a, mid + 1, hi);
	if (thread_a[mid] > thread_a[mid + 1]) {//if already not sorted, merge subarrays

		merge(thread_a, thread_tmp, lo, mid, hi);
	}
	else {
		unsigned int i;
		for (i = lo; i <= hi; i++)
			thread_tmp[i] = thread_a[i];
	}



}

void* thread_merge (void* my_rank){ 
	/*4 subarrays already sorted in "tmp",
		2 larger sorted subarrays are made in "a"
	*/
	unsigned int rank = (unsigned int) my_rank;
	unsigned int lo = rank*(n/2);
	unsigned int hi = (rank+1)*(n/2) - 1;
	unsigned int mid = lo + (hi-lo)/2 ;
	
	unsigned int i, k, j;
	/*
	for(k=lo; k<=hi; k++){
		tmp[k] = a[k];
	}
	*/
	i = lo;
	j = mid+1;
	for(k=lo; k<=hi; k++){
		if(i>mid){
			a[k] = tmp[j];
			j++;
		}else if(j>hi){
			a[k] = tmp[i];
			i++;
		}else if(tmp[j] < tmp[i]){
			a[k] = tmp[j];
			j++;
		}else{
			a[k] = tmp[i];
			i++;
		}
	}
	return NULL;
	
}

void* thread_sort(void* my_rank){
	unsigned int rank = (unsigned int) my_rank;
	unsigned int lo = rank*(n/4);
	unsigned int hi = (rank+1)*(n/4) - 1;
	unsigned int mid = lo + (hi-lo)/2 ;
	
	unsigned int i;
	
	for(i=lo; i<=hi;i++){ // copy "a" in "tmp" by 4 threads in parallel
		tmp[i] = a[i];
	}
	
	if(hi <= lo) return NULL;
	mergeSort(tmp, a, lo, mid);
	mergeSort(tmp, a, mid + 1, hi);
	if (a[mid] <= a[mid + 1]) { //if already sorted
		for (i = lo; i <= hi; i++)
			tmp[i] = a[i];
		return NULL;
	}
	
	merge(a, tmp, lo, mid, hi);
	return NULL;
}

void* thread_search(void* input){
	struct Search* s = (struct Search*) input;
	
	unsigned int lo = (*s).lo;
	unsigned int hi = (*s).hi;
	unsigned int rank = (*s).rank;
	int num = (*s).num;
	

	/*
	if(rank<2){//find position in second half
		lo = (n/2);
		hi = n - 1;	
	}else{ // find position in first half
		lo = 0;
		hi = (n/2) - 1;
	}
	
	int num = a[rank*(n/4)];
	*/
	
	if(num < a[lo]){
		samples[rank] = lo;
		return NULL;
	}
	if(num > a[hi]){
		samples[rank] = hi;
		return NULL;
	}
	
	while (lo <= hi) {
		unsigned int mid = lo + (hi-lo)/2;
 
        if(a[mid] == num){
            samples[rank] = mid;
			return NULL;
		}else if(a[mid] < num){
            lo = mid + 1;
		}else{
            hi = mid - 1;
		}
	}
    // insert position
	samples[rank] = hi+1;
    return NULL;
}

void* thread_mergeParallel(void* input){
	struct Bound* b = (struct Bound*) input;
	unsigned int lo1 =(*b).lo1;
	unsigned int hi1 =(*b).hi1;
	unsigned int lo2 =(*b).lo2;
	unsigned int hi2 =(*b).hi2;
	unsigned int offset =(*b).offset;
	unsigned int length =(*b).length;
	
	unsigned int i,j,k;
	
	i = lo1;
	j = lo2;
	for(k=offset; k< length+offset; k++){
		if(i>=hi1){
			tmp[k] = a[j];
			j++;
		}else if(j>=hi2){
			tmp[k] = a[i];
			i++;
		}else if(a[j] < a[i]){
			tmp[k] = a[j];
			j++;
		}else{
			tmp[k] = a[i];
			i++;
		}
	}
	
	return NULL;
	
}	

void mergeSortParallel (const int* values, unsigned int N, int* sorted) {
	//tmp=(int*) malloc(N * sizeof(int));
	tmp = sorted;
    n=N;
    a=(int*)values;
	
	pthread_t* handles = (pthread_t*) malloc (4*sizeof(pthread_t));

	unsigned int j;
	//stage 1 : sort 4 subarrays of "a" in "tmp" using merge sort
	
	for (j=0; j<4; j++){  
      pthread_create( &handles[j], NULL, thread_sort, (void*)j );  
	}
	for (j=0; j<4; j++) {
      pthread_join( handles[j], NULL); 
	}
	
 
	//stage 2 : use 2 threads to merge 4 subarrays of "tmp" into 2 larger subarrays in "a"
	
	for (j=0; j<2; j++){  
      pthread_create( &handles[j], NULL, thread_merge, (void*)j );  
	}
	for (j=0; j<2; j++) {
      pthread_join( handles[j], NULL); 
	}

	//stage 3 : search samples in the other half
	//////////////////////////////////////////////////
	struct Search s[4];
	for(j=0;j<4;j++){//limit search length from n/2 to n/4
		if(j<2){//thread 0 , 1 find index in second half
			if(a[j*(n/4)]<a[3*(n/4)]){
				s[j].lo = n/2;
				s[j].hi = 3*(n/4) - 1;
			}else{
				s[j].lo = 3*(n/4);
				s[j].hi = n-1;
			}
			
		}else{//thread 2 , 3 find index in first half
			if(a[j*(n/4)]<a[(n/4)]){
				s[j].lo = 0;
				s[j].hi = (n/4) - 1;
			}else{
				s[j].lo = n/4;
				s[j].hi = (n/2) - 1;
			}
		}
		s[j].rank = j;
		s[j].num = a[j*(n/4)];
	}
	//////////////////////////////////////////////////
	
	for (j=0; j<4; j++){  
      pthread_create( &handles[j], NULL, thread_search, (void*)(&s[j]) );  
	}
	for (j=0; j<4; j++) {
      pthread_join( handles[j], NULL); 
	}

	//stage 4 : sort boundaries in each half
	//sort 0,samples[2],n/4,samples[3],n/2
	
	boudary1[0] = 0;	
	boudary1[4] = n/2;	
	if(samples[2] > n/4){
		boudary1[1] = n/4;
		boudary1[2] = samples[2];
	}else{
		boudary1[1] = samples[2];
		boudary1[2] = n/4;
	}
	if(boudary1[2]>samples[3]){
		boudary1[3] = boudary1[2];
		boudary1[2] = samples[3];
	}else{
		boudary1[3] = samples[3];
	}
	
	//sort n/2,samples[0],3n/4,samples[1],n
	
	boudary2[0] = n/2;
	boudary2[4] = n;
	if(samples[0] > 3*(n/4)){
		boudary2[1] = 3*(n/4);
		boudary2[2] = samples[0];
	}else{
		boudary2[1] = samples[0];
		boudary2[2] = 3*(n/4);
	}
	if(boudary2[2]>samples[1]){
		boudary2[3] = boudary2[2];
		boudary2[2] = samples[1];
	}else{
		boudary2[3] = samples[1];
	}
	
	
	//stage 5 : merge in parallel from "a" to "tmp"
	
	//set boundaries in first and second half for each thread
	
	struct Bound bounds [4];
	for(j=0; j<4; j++){
		bounds[j].lo1 = boudary1[j];
		bounds[j].hi1 = boudary1[j+1];
		bounds[j].lo2 = boudary2[j];
		bounds[j].hi2 = boudary2[j+1];
		bounds[j].offset = boudary1[j]+boudary2[j]-(n/2);
		bounds[j].length = (boudary1[j+1]-boudary1[j])+(boudary2[j+1]-boudary2[j]);
	}
	
	
	for (j=0; j<4; j++){  
      pthread_create( &handles[j], NULL, thread_mergeParallel, (void*)(&bounds[j]) );  
	}
	for (j=0; j<4; j++) {
      pthread_join( handles[j], NULL); 
	}

	//copy sorted values to "sorted"
	//no copy needed since array "tmp" is array "sorted" itself
	/*
	for(j=0; j<n; j++){
		sorted[j] = tmp[j];
	}
	
	free(tmp);
	*/
	
	
	
	
}