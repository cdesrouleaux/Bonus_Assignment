//
// Jean Desrouleaux
// Programming Assign 3 - movieline.c
// 03/09/24
// Dr. Neslisah Torosdagli
// COP3505 Tu Th 6-7:15pm
//

#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

int extraMemoryAllocated;

void *Alloc(size_t sz)
{
	extraMemoryAllocated += sz;
	size_t* ret = malloc(sizeof(size_t) + sz);
	*ret = sz;
	printf("Extra memory allocated, size: %ld\n", sz);
	return &ret[1];
}

void DeAlloc(void* ptr)
{
	size_t* pSz = (size_t*)ptr - 1;
	extraMemoryAllocated -= *pSz;
	printf("Extra memory deallocated, size: %ld\n", *pSz);
	free((size_t*)ptr - 1);
}

size_t Size(void* ptr)
{
	return ((size_t*)ptr)[-1];
}

void swap(int* x, int* y)
{
	int temp = *x;
	*x = *y;
	*y = temp;
}

// O(n) runtime faster than how you typically create a max heap which is O(nlogn)
void f_heap(int *pData, int n, int i)
{
    // make sure to make 2 if's int 1 if 1 else to to avoid cheating
    int lg=i;
    int l = 2*i, r = 2*i+1;
    // executes if total tree size is bigger that left child
    if(n > l)
    {
        // executes if parent is smaller than left child and makes child the new parent
        if(pData[lg] < pData[l])
            lg = l; 
    } 
    // executes if total tree size is bigger that right child
    if(n > r)
    {   
        // executes if parent is smaller than right child and makes child the new parent
        if(pData[lg] < pData[r])
            lg = r;
    }

    if(i != lg)
    {
        swap(&pData[lg], &pData[i]);
        f_heap(pData, n, lg);  
    }      
	
}

// implements heap sort
// extraMemoryAllocated counts bytes of memory allocated
// heap sorting uses the heap/ tree structure to sort. First you must obtain the max heap then you can sort.
void heapSort(int pData[], int n)
{
    int s = 0;
	// for loop first creates a max heap
	for (int i = floor(n/2); i > -1; i--)
		f_heap(pData, n, i);

	// then the following for loop takes the max heap and sorts it
   for (int i = n-1; i > -1; i--)
   {
        swap(&pData[s], &pData[i]);
        f_heap(pData, i, s);
   }
}

void merge(int pData[], int l, int mid, int r)
{
    int z1 = mid - l + 1;
    int z2 =  r - mid;

    int *L = (int*) malloc(z1*sizeof(int));
    int *R = (int*) malloc(z2*sizeof(int));

    for (int h1 = 0; h1 < z1; h1++)
        L[h1] = pData[l + h1];
    for (int h2 = 0; h2 < z2; h2++)
        R[h2] = pData[mid + 1 + h2];


    int i=0, j=0, k=l;

    while (i < z1 && j < z2)
    {
        if (L[i] <= R[j])
        {
            pData[k] = L[i];
            i++;
        }
        else
        {
            pData[k] = R[j];
            j++;
        }
        k++;
    }

    while (i < z1)
    {
        pData[k] = L[i];
        i++;
        k++;
    }

    while (j < z2)
    {
        pData[k] = R[j];
        j++;
        k++;
    }
    free(L);
    free(R);
}

// merge sort splits data into halves until you get singles than data is reassembled into a sorted array
// implement merge sort
// extraMemoryAllocated counts bytes of extra memory allocated
void mergeSort(int pData[], int l, int r)
{
    if(l < r)
    {
        int mid = (l+r) / 2; // middle address of pData array

        mergeSort(pData, l, mid); // recursivly halves the left half of the array
        mergeSort(pData, mid+1, r); // recursivly halves the right half of the array

        merge(pData, l, mid, r); // merges both halves/ subarrays of pData
    }
}

// implement insertion sort
// extraMemoryAllocated counts bytes of memory allocated
// insertion pulls the data out from its current slot and places it to the left into its correct slot and slides over the remaining data to the right 1
void insertionSort(int* pData, int n)
{
    int item,j;
    for(int i = 1; i < n; i++)
    {
        item = pData[i];
        for(j = i-1; j >= 0; j--)
        {
            if(pData[j] > item)
                pData[j+1] = pData[j];

            else
                break;
        }
        // j always ends up as the smallest in first position or on the right of first position
        pData[j+1] = item;
    }
}

// implement bubble sort
// extraMemoryAllocated counts bytes of extra memory allocated
// bubble sort swaps adjacent elements repeatedly until entire array is sorted takes the longest
void bubbleSort(int* pData, int n)
{
    int flag;

    //swaps adjacent values continuously until the entire array is sorted
    while(1)
    {    
        flag = 0;

        for(int i = 0; i < n; i++)
        {
            // if left is bigger than right swap
            if(pData[i+1] < pData[i])
            {
                swap(&pData[i], &pData[i+1]);
                flag = 1;
            }
        }
            /*breaks if none of the values need to be swapped (flag=0). This means that the 
            left values are always smaller than the right values and the array is sorted*/
        if(flag == 0) break;
    }
}

// implement selection sort
// extraMemoryAllocated counts bytes of extra memory allocated
// finds the smallest element in the array is sticks it in front of the search space repeatedly until entire array is sorted
void selectionSort(int* pData, int n)
{
	int  min_index;
	// greater for loop shrinks the search space [0,..,10] [1,..,10] [2,..,10] [9,..,10]
	for(int i=0; i < n-1; i++)
	{
		//sets min to the first value of the search space
		min_index = i;
		for(int j=i+1; j < n; j++) //sets j to the adjacent value in the search space
		{
			//finds the smallest value
			if(pData[j] < pData[min_index])
				min_index = j;		 
		}
		//puts the smallest value at the front of the search space
            swap(&pData[i], &pData[min_index]);
	}
}

// parses input file to an integer array
int parseData(char *inputFileName, int **ppData)
{
	FILE* inFile = fopen(inputFileName,"r");
	int dataSz = 0;
	*ppData = NULL;
	
	if (inFile)
	{
		fscanf(inFile,"%d\n",&dataSz);
		*ppData = (int *)Alloc(sizeof(int) * dataSz);
		// Implement parse data block
	}
	
	return dataSz;
}

// prints first and last 100 items in the data array
void printArray(int pData[], int dataSz)
{
	int i, sz = dataSz - 100;
	printf("\tData:\n\t");
	for (i=0;i<100;++i)
	{
		printf("%d ",pData[i]);
	}
	printf("\n\t");
	
	for (i=sz;i<dataSz;++i)
	{
		printf("%d ",pData[i]);
	}
	printf("\n\n");
}

int main(void)
{
	clock_t start, end;
	int i;
    double cpu_time_used;
	char* fileNames[] = {"input1.txt", "input2.txt", "input3.txt"};
	
	for (i=0;i<3;++i)
	{
		int *pDataSrc, *pDataCopy;
		int dataSz = parseData(fileNames[i], &pDataSrc);
		
		if (dataSz <= 0)
			continue;
		
		pDataCopy = (int *)Alloc(sizeof(int)*dataSz);
	
		printf("---------------------------\n");
		printf("Dataset Size : %d\n",dataSz);
		printf("---------------------------\n");
		
		printf("Selection Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		selectionSort(pDataCopy, dataSz);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);
		
		printf("Insertion Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		insertionSort(pDataCopy, dataSz);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);

		printf("Bubble Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		bubbleSort(pDataCopy, dataSz);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);
		
		printf("Merge Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		mergeSort(pDataCopy, 0, dataSz - 1);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);

        printf("Heap Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		heapSort(pDataCopy, dataSz - 1);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);
		
		DeAlloc(pDataCopy);
		DeAlloc(pDataSrc);
	}
	
}
