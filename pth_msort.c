// Include your C header files here
#include "pth_msort.h"
#include <stdlib.h>
#include <stdio.h>

unsigned int size;

// define a struct to package sort module arguments 
struct SortArgs 
{
	int* data;
	unsigned int lower;
	unsigned int upper;
};

// define a struct to package merge module arguments 
struct MergeArgs 
{
	int* data;
	unsigned int lower;
	unsigned int upper;
};

// define a struct to package parallel merge module arguments 
struct PMergeArgs 
{
	int* data;
	int* sorted;
	unsigned int* partitions;
	unsigned char part;
};

// swap function used in quickSort algorithm
void swap(int* a, int* b)
{
    int t = *a;
    *a = *b;
    *b = t;
}

// partition function used in quickSort algorithm
int partition(int* data, int low, int high)
{
    int pivot = data[high]; 
    int i = (low- 1); 
    int j;
	
    for (j = low; j <= high - 1; j++) {
        if (data[j] < pivot) {
            i++; 
            swap(&data[i], &data[j]);
        }
    }
    swap(&data[i + 1], &data[high]);
    return (i + 1);
}

// quickSort algorithm
void quickSort(int* data, int low, int high)
{
    if (low < high) 
	{
        int pi = partition(data, low, high);
        quickSort(data, low, pi - 1);
        quickSort(data, pi + 1, high);
    }
}

// sort module
void* SortModule(void* arguments)
{
	struct SortArgs* args = (struct SortArgs*) arguments;	
	quickSort(args->data, args->lower, args->upper);             	
}

// merge algorithm
void Merge(int* data, unsigned int lower, unsigned int middle, unsigned int upper)
{	
	// define loop counters
	unsigned int i, j, k;
	// determine the size of data parts
	int n1 = middle - lower + 1;
    int n2 = upper - middle;
	// define a variable to keep the lower part of data
	int* l_data = (int*) malloc(n1*sizeof(int));
	// define a variable to keep the upper part of data
	int* u_data = (int*) malloc(n2*sizeof(int));	
		
	// move the lower part of the data to l_data
	for(i = 0; i < n1; i++)
		l_data[i] = data[lower + i];
	// move the upper part of the data to u_data			
	for(j = 0; j < n2; j++)
		u_data[j] = data[middle + 1 + j];
	// main merge algorithm
	i = 0;
	j = 0;
	k = lower;
	
	while (i < n1 && j < n2) 
    {
        if (l_data[i] <= u_data[j]) 
            data[k++] = l_data[i++];
        else      
            data[k++] = u_data[j++];
    }
	
	while (i < n1) 
        data[k++] = l_data[i++];
    
    while (j < n2)   
        data[k++] = u_data[j++];
    
	// set allocated space free
	free(l_data);
	free(u_data);
}

// merge module
void* MergeModule(void* arguments)
{
	struct MergeArgs* args = (struct MergeArgs*) arguments;	
	unsigned int middle = args->lower + (args->upper - args->lower)/2;
	Merge(args->data, args->lower, middle, args->upper);             	
}

// binary search algorithms 
unsigned int binarySearchCount(int* data, unsigned int upper, int x)
{
	unsigned int count = 0, mid = 0, lower = 0;
	while (lower <= upper) 
	{
		mid = (upper + lower) / 2;
		if (data[mid] <= x) 
		{
			count = mid + 1;
			lower = mid + 1;
		}
		else
		upper = mid - 1;
	}
	return count;
}

// find partitions of data for parallel merge
void FindPartitions(int* data, unsigned int* partitions)
{
	// find the array length 
	unsigned int M = size/2;
	// define a loop counter
	unsigned int i;
	
	// split the data

	int* u_data = (int*) malloc(M*sizeof(int));	
		
	// move the lower part of the data to l_data
	for(i = 0; i < M; i++)	
		u_data[i] = data[M + i];
	
	
	// find partitions
	partitions[0] = 0;
	partitions[1] = binarySearchCount(u_data, M-1, data[(M/4) - 1]);
	partitions[2] = binarySearchCount(u_data, M-1, data[(M/2) - 1]);
	partitions[3] = binarySearchCount(u_data, M-1, data[(3*M/4) - 1]);
	partitions[4] = M;
	
	// set allocated space free
	free(u_data);	
}
 
// parallel merge algorithm
void PMerge(int* data, int* sorted, unsigned int* partitions, unsigned char part)
{	
	// define loop counters
	unsigned int i, j, k;
	// determine array sizes
	int n1 = size/8;
    int n2 = partitions[part+1] - partitions[part];
	
	// move the lower part of the data to l_data
	int sp1 = part*n1;
	
	// main merge algorithm
	if ( n2 == 0 )
	{
		for(i = 0; i < n1; i++)
			sorted[sp1 + partitions[part] + i] = data[sp1 + i];
	}
	else
	{
		i = 0;
		j = 0;
		k = 0;
		while (i < n1 && j < n2) 
		{
			if (data[sp1+i] <= data[size/2+partitions[part]+j]) 
			{
				sorted[sp1 + partitions[part]+k] = data[sp1+i];
				i++;
			}				
			else  
			{				
				sorted[sp1 + partitions[part]+k] = data[size/2+partitions[part]+j];
				j++;			
			}
			k++;
		}
		
		while (i < n1)
		{		
			sorted[sp1 + partitions[part]+k] = data[sp1+i];
			i++;
			k++;
		}
		
		while (j < n2) 
		{		
			sorted[sp1 + partitions[part]+k] = data[size/2+partitions[part]+j];
			j++;
			k++;
		}
	}	
}

// parallel merge module
void* PMergeModule(void* arguments)
{
	struct PMergeArgs* args = (struct PMergeArgs*) arguments;	
	PMerge(args->data, args->sorted, args->partitions, args->part);             	
}

// main code
void mergeSortParallel(const int* values, unsigned int N, int* sorted) 
{	
	// define variables to keep data parts
	int* data = (int*) values;
	// define four thread handles
	pthread_t* handles = (pthread_t*) malloc (4*sizeof(pthread_t));
	// define a loop counter
	unsigned char i;
	// define a variable to keep partitions
	unsigned int* partitions = (unsigned int*) malloc (5*sizeof(unsigned int));; 
	
	// fill global variables
	size = N;
	unsigned int hsize = N/2;
	unsigned int qsize = N/4;


	// call 4 threads where each thread calls a sort module		
	struct SortArgs arguments1 = {data, 0, qsize - 1};
	struct SortArgs arguments2 = {data, qsize, hsize - 1};
	struct SortArgs arguments3 = {data, hsize, 3*qsize - 1};
	struct SortArgs arguments4 = {data, 3*qsize, N - 1};	
	
	pthread_create(&handles[0], NULL, SortModule, (void*)(&arguments1));
	pthread_create(&handles[1], NULL, SortModule, (void*)(&arguments2));
	pthread_create(&handles[2], NULL, SortModule, (void*)(&arguments3));
	pthread_create(&handles[3], NULL, SortModule, (void*)(&arguments4));
	
	// wait for called threads to finish
	for (i=0; i<4; i++) 
		pthread_join(handles[i], NULL);
		
	// call 2 threads where each thread calls a merge module
	struct MergeArgs arguments5 = {data, 0, hsize - 1};
	struct MergeArgs arguments6 = {data, hsize, N - 1};
	
	pthread_create(&handles[0], NULL, MergeModule, (void*)(&arguments5));
	pthread_create(&handles[1], NULL, MergeModule, (void*)(&arguments6));
	
	// wait for called threads to finish
	for (i=0; i<2; i++) 
		pthread_join(handles[i], NULL);
	
	// find partitions of parallel merge algorithm
	FindPartitions(data, partitions);
	
	// call 4 threads to do the final merge
	struct PMergeArgs pArg0 = {data, sorted, partitions, 0};
	struct PMergeArgs pArg1 = {data, sorted, partitions, 1};
	struct PMergeArgs pArg2 = {data, sorted, partitions, 2};
	struct PMergeArgs pArg3 = {data, sorted, partitions, 3};

	pthread_create(&handles[0], NULL, PMergeModule, (void*)(&pArg0));
	pthread_create(&handles[1], NULL, PMergeModule, (void*)(&pArg1));
	pthread_create(&handles[2], NULL, PMergeModule, (void*)(&pArg2));
	pthread_create(&handles[3], NULL, PMergeModule, (void*)(&pArg3));
	
	// wait for called threads to finish
	for (i=0; i<4; i++) 
		pthread_join(handles[i], NULL);
			
	// set allocated space free
	free(partitions);
	// set thread handles free
	free(handles);	
}
