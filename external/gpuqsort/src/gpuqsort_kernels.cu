/**
* GPU Quicksort Kernels
* ---------------------
* Functions part1, part2, part3 and lqsort
* Copyright 2007-2008 Daniel Cederman and Philippas Tsigas
*
* This work is licensed under the Creative Commons
* Attribution-Noncommercial-No Derivative Works 3.0
* Unported License. To view a copy of this license,
* visit http://creativecommons.org/licenses/by-nc-nd/3.0/
* or send a letter to Creative Commons, 171 Second Street,
* Suite 300, San Francisco, California, 94105, USA.
*
* ----
*
* Functions bitonicSort and cumcount are based upon the
* NVIDIA CUDA SDK Code Samples "Bitonic Sort" and "Scan"
*
* "This code is released free of charge for use in
* derivative works, whether academic, commercial, or personal"
* - CUDA Website
* 
* http://developer.download.nvidia.com/licenses/general_license.txt
*
**/

#include "gpuqsort.h"

#undef THREADS
#define THREADS blockDim.x

extern __shared__ unsigned int sarray[]; 

#ifdef HASATOMICS
__device__ unsigned int ohtotal = 0;
#endif

/**
* Swaps the location of two unsigned ints
* @param a This unsigned int will swap place with unsigned int b
* @param b This unsigned int will swap place with unsigned int a
*/
//template <typename unsigned int>
__device__ inline void swap(unsigned int& a, unsigned int& b)
{
    unsigned int tmp = a;
    a = b;
    b = tmp;
}

/**
* Perform a bitonic sort
* @param values The unsigned ints to be sorted
* @param target Where to place the sorted unsigned int when done
* @param size The number of unsigned ints
*/
//template <typename unsigned int>
__device__ inline
void bitonicSort(unsigned int* fromvalues, unsigned int* tovalues, unsigned int from, unsigned int size)
{
	unsigned int* shared = (unsigned int*)sarray;

	unsigned int coal = (from&0xf);
	size = size + coal;
	from = from - coal;

	int sb = 2 << (int)(__log2f(size));

	// Buffer data to be sorted in the shared memory
	for(int i=threadIdx.x;i<size;i+=THREADS)
	{
		shared[i] = fromvalues[i+from];
	}

	for(int i=threadIdx.x;i<coal;i+=THREADS)
		shared[i]=0;

	// Pad the data
	for(int i=threadIdx.x+size;i<sb;i+=THREADS)
		shared[i] = 0xffffffff;

    __syncthreads();

    // Parallel bitonic sort.
    for (int k = 2; k <= sb; k *= 2)
    {
        // Bitonic merge:
        for (int j = k / 2; j>0; j /= 2)
        {
			for(int tid=threadIdx.x;tid<sb;tid+=THREADS)
			{
				unsigned int ixj = tid ^ j;
	            
				if (ixj > tid)
				{
					if ((tid & k) == 0)
					{
						if (shared[tid] > shared[ixj])
						{
							swap(shared[tid], shared[ixj]);
						}
					}
					else
					{
						if (shared[tid] < shared[ixj])
						{
							swap(shared[tid], shared[ixj]);
						}
					}
				}
            }
            
            __syncthreads();
        }
    }
	__syncthreads();

	// Write back the sorted data to its correct position
	for(int i=threadIdx.x;i<size;i+=THREADS)
		if(i>=coal)
		    tovalues[i+from] = shared[i];
	__syncthreads();
}


/**
* Perform a cumulative count on two arrays
* @param lblock Array one
* @param rblock Array two
*/
__device__ inline void cumcount(unsigned int *lblock, unsigned int *rblock)
{
	int tx = threadIdx.x;

    int offset = 1; 
 
    __syncthreads();

	for (int d = THREADS>>1; d > 0; d >>= 1) // build sum in place up the tree 
    {
        __syncthreads();

		if (tx < d)    
        { 
			int ai = offset*(2*tx+1)-1;
            int bi = offset*(2*tx+2)-1;
            lblock[bi] += lblock[ai];
			rblock[bi] += rblock[ai];
		} 
        offset *= 2; 
    } 
	__syncthreads(); 
    if (tx == 0) 
	{ 
		lblock[THREADS] = lblock[THREADS-1];
		rblock[THREADS] = rblock[THREADS-1];
		lblock[THREADS - 1] =0;
		rblock[THREADS - 1] =0; 
	} // clear the last unsigned int */
	__syncthreads(); 

    for (int d = 1; d < THREADS; d *= 2) // traverse down tree & build scan 
    { 
        offset >>= 1; 
        __syncthreads(); 
 
        if (tx < d) 
        { 
			int ai = offset*(2*tx+1)-1; 
            int bi = offset*(2*tx+2)-1; 
 
            int t   = lblock[ai]; 
			lblock[ai]  = lblock[bi]; 
            lblock[bi] += t; 

            t   = rblock[ai]; 
            rblock[ai]  = rblock[bi]; 
            rblock[bi] += t; 

        } 
    } 
}


/**
* Part One - Counts the number of unsigned ints larger or smaller than the pivot. It then
* performs a cumulative sum so that each thread knows where to write
* @param data   unsigned ints to be counted
* @param params Specifies which data each thread block is responsible for
* @param hist   The cumulative sum for each thread is stored here
* @param lengths The total sum for each thread block is stored here
*/
//template <typename unsigned int>
__global__ void part1(unsigned int* data, struct Params<unsigned int>* params, struct Hist* hist, Length<unsigned int>* lengths)
{
	const int tx = threadIdx.x;

	unsigned int* lblock = (unsigned int*)sarray;
	unsigned int* rblock = (unsigned int*)(&lblock[(blockDim.x+1)]);
	unsigned int* minpiv = (unsigned int*)(&rblock[(blockDim.x+1)]);
	unsigned int* maxpiv = (unsigned int*)(&minpiv[blockDim.x]); 


	// Where should we read?
	unsigned int start = params[blockIdx.x].from;
	unsigned int end = params[blockIdx.x].end;
	unsigned int pivot = params[blockIdx.x].pivot;

	// Stores the max and min value of the data. Used to decide a new pivot
	minpiv[tx] = data[start+tx];
	maxpiv[tx] = data[start+tx];

	__syncthreads();
	int ll=0;
	int lr=0;

	__syncthreads();


	int coal = (start&0xf);
	start = start-coal;

	// Go through the data
	if(tx+start<end)
	{
		unsigned int d = data[tx+start];

		if(!(tx<coal))
		{

		// Counting unsigned ints smaller...
		if(d<pivot)
			ll++;
		else
		// or larger than the pivot
		if(d>pivot)
			lr++;

		// Store the max and min unsigned int
		minpiv[tx] = min(minpiv[tx],d);
		maxpiv[tx] = max(maxpiv[tx],d);
		}
	}


	// Go through the data
	for(unsigned int i=tx+start+THREADS;i<end;i+=THREADS)
	{
		unsigned int d = data[i];

		// Counting unsigned ints smaller...
		if(d<pivot)
			ll++;
		else
		// or larger than the pivot
		if(d>pivot)
			lr++;

		// Store the max and min unsigned int
		minpiv[tx] = min(minpiv[tx],d);
		maxpiv[tx] = max(maxpiv[tx],d);
	}

	lblock[tx]=ll;
	rblock[tx]=lr;

	__syncthreads();

	// Perform a cumulative sum
    cumcount((unsigned int*)lblock,(unsigned int*)rblock);

    if(tx==0)
    {
		// Decide on max and min unsigned int
		for(int i=0;i<THREADS;i++)
		{
			minpiv[0]=min(minpiv[0],minpiv[i]);
			maxpiv[0]=max(maxpiv[0],maxpiv[i]);
		}
	  }
  	__syncthreads();
	
	// Store each threads part of the cumulative count
	hist->left[blockIdx.x*(THREADS)+threadIdx.x]  = lblock[threadIdx.x+1];
	hist->right[blockIdx.x*(THREADS)+threadIdx.x] = rblock[threadIdx.x+1];

	// Store the total sum
	lengths->left[blockIdx.x]  = lblock[THREADS];
	lengths->right[blockIdx.x] = rblock[THREADS];

	// Store the max and min unsigned int
	lengths->minpiv[blockIdx.x] = minpiv[0];
	lengths->maxpiv[blockIdx.x] = maxpiv[0];

}

/**
* Part Two - Move unsigned ints to their correct position in the auxillary array
* @param data   unsigned ints to be moved
* @param data2  Destination for unsigned ints
* @param params Specifies which data each thread block is responsible for
* @param hist   The cumulative sum for each thread is stored here
* @param lengths The total sum for each thread block is stored here
*/
//template <typename unsigned int>
__global__ void part2(unsigned int* data, unsigned int* data2, struct Params<unsigned int>* params, struct Hist* hist, Length<unsigned int>* lengths)
{
	const int tx = threadIdx.x;
	const int bx = blockIdx.x;

	// Each thread uses the cumulative sum to know where to write
	unsigned int x = lengths->left[bx] + hist->left[bx*(THREADS)+tx]-1;// - 1;
	unsigned int y = lengths->right[bx] - hist->right[bx*(THREADS)+tx];

	// Where should we read?
	unsigned int start = params[bx].from;
	unsigned int end = params[bx].end;
	unsigned int pivot = params[bx].pivot;

	__syncthreads();

	int coal = (start&0xf);
	start = start-coal;

	// Go through all the assigned data
	if(tx+start<end)
	{
		// Reading unsigned ints...
		unsigned int d = data[tx+start];

		if(!(tx<coal))
		{

		// and writing them to auxillary array
		if(d<pivot)
			data2[x--]=d;
		else
		if(d>pivot)
			data2[y++]=d;
		}
	}

	__syncthreads();

	// Go through all the assigned data
	for(unsigned int i=start+tx+THREADS;i<end;i+=THREADS)
	{
		// Reading unsigned ints...
		unsigned int d = data[i];

		// and writing them to auxillary array
		if(d<pivot)
		{
			data2[x--]=d;
		}
		else
		if(d>pivot)
			data2[y++]=d;
	}

	return;
}

/**
* Part Three - Write the pivot value
* @param data   Destination for pivot
* @param params Specifies which data each thread block is responsible for
* @param hist   The cumulative sum for each thread is stored here
* @param lengths The total sum for each thread block is stored here
*/
//template <typename unsigned int>
__global__ void part3(unsigned int* data, struct Params<unsigned int>* params, struct Hist* hist, Length<unsigned int>* lengths)
{
	const int tx = threadIdx.x;
	const int bx = blockIdx.x;

	// If we are the "last" thread block that is assigned to the same data sequence
	// we write the pivot between the left and right block
	if(params[bx].last)
	{
		// Get destination position
		unsigned int x = lengths->left[bx] + hist->left[bx*THREADS+THREADS-1] + tx;
		unsigned int y = lengths->right[bx] - hist->right[bx*THREADS+THREADS-1];
		unsigned int pivot = params[bx].pivot;

		// Write the pivot values
		for(;x<y;x+=THREADS)
			data[x]=pivot;
	}
}

/**
* The local quicksort - sorts a block of data with no inter-block synchronization
* @param adata  Contains some of the blocks to be sorted and also acts as the final
*               destination for sorted data
* @param adata2 Contains some of the blocks to be sorted
* @param bs     List of blocks to be sorted and a pointer telling if a specific block is
*               in \a adata or \a adata2
*/
//template <typename unsigned int>
__global__ void lqsort(unsigned int* adata, unsigned int* adata2, struct LQSortParams* bs, unsigned int phase)
{
	__shared__ unsigned int lphase;
	lphase=phase;

	// Shorthand for the threadid
	int tx = threadIdx.x;

	// Stack pointer
	__shared__ int bi;
	
	// Stack unsigned ints
	__shared__ unsigned int beg[32];
	__shared__ unsigned int end[32];
	__shared__ bool flip[32];

	unsigned int* lblock = (unsigned int*)sarray;
	unsigned int* rblock = (unsigned int*)(&lblock[(blockDim.x+1)]);


	// The current pivot
	__shared__ unsigned int pivot;

	// The sequence to be sorted
	__shared__ unsigned int from;
	__shared__ unsigned int to;

	// Since we switch between the primary and the auxillary buffer,
	// these variables are required to keep track on which role
	// a buffer currently has
	__shared__ unsigned int* data;
	__shared__ unsigned int* data2;
	__shared__ unsigned int sbsize;

	__shared__ unsigned int bx;
	if(threadIdx.x==0)
#ifdef HASATOMICS
		bx = atomicInc(&ohtotal,50000);
#else
		bx = blockIdx.x;
#endif
		
	__syncthreads();

	while(bx<gridDim.x)
	{

	

	// Thread 0 is in charge of the stack operations
	if(tx==0)
	{
		// We push our first block on the stack
		// This is the block given by the bs parameter
		beg[0] = bs[bx].beg;
		end[0] = bs[bx].end;
		flip[0] = bs[bx].flip;
		sbsize = bs[bx].sbsize;

		bi = 0;
	}

	__syncthreads();

	// If we were given an empty block there is no need to continue
	if(end[0]==beg[0])
		return;

	// While there are items left on the stack to sort
	while(bi>=0)
	{
		__syncthreads();
		// Thread 0 pops a fresh sequence from the stack
		if(tx==0)
		{
			from = beg[bi];
			to = end[bi];

			// Check which buffer the sequence is in
			if(!flip[bi])
			{
				data = adata2;
				data2 = adata;
			}
			else
			{
				data = adata;
				data2 = adata2;
			}

		}
	

		__syncthreads();

		// If the sequence is smaller than SBSIZE we sort it using
		// an alternative sort. Otherwise each thread would sort just one
		// or two unsigned ints and that wouldn't be efficient
		if((to-from)<(sbsize-16))
		{
			// Sort it using bitonic sort. This could be changed to some other
			// sorting method. Store the result in the final destination buffer
			if((to-from>=1)&&(lphase!=2))
				bitonicSort(data,adata,from,to-from);
			__syncthreads();

			// Decrement the stack pointer
			if(tx==0)
				bi--;
			__syncthreads();
			// and continue with the next sequence
			continue;
		}

 
		if(tx==0)
		{
			// Create a new pivot for the sequence
			// Try to optimize this for your input distribution
			// if you have some information about it
			unsigned int mip = min(min(data[from],data[to-1]),data[(from+to)/2]);
			unsigned int map = max(max(data[from],data[to-1]),data[(from+to)/2]);
			pivot = min(max(mip/2+map/2,mip),map);
		}


		unsigned int ll=0;
		unsigned int lr=0;

		__syncthreads();
		
		unsigned int coal = (from)&0xf;

		if(tx+from-coal<to)
		{
			unsigned int d = data[tx+from-coal];

			if(!(tx<coal))
			{
				// Counting unsigned ints that have a higher value than the pivot
				if(d<pivot)
					ll++;
				else
				// or a lower
				if(d>pivot)
					lr++;
			}
		}


		// Go through the current sequence
		for(int i=from+tx+THREADS-coal;i<to;i+=THREADS)
		{
			unsigned int d = data[i];

			// Counting unsigned ints that have a higher value than the pivot
			if(d<pivot)
				ll++;
			else
			// or a lower
			if(d>pivot)
				lr++;
		}

		// Store the result in a shared array so that we can calculate a
		// cumulative sum
		lblock[tx]=ll;
		rblock[tx]=lr;
		
		__syncthreads();

		// Calculate the cumulative sum
		cumcount((unsigned int*)lblock,(unsigned int*)rblock);

		__syncthreads();

		// Let thread 0 add the new resulting subsequences to the stack
		if(tx==0)
		{
			// The sequences are in the other buffer now
			flip[bi+1] = !flip[bi];
			flip[bi] = !flip[bi];

			// We need to place the smallest object on top of the stack
			// to ensure that we don't run out of stack space
			if(lblock[THREADS]<rblock[THREADS])
			{
				beg[bi+1]=beg[bi];
				beg[bi]=to-rblock[THREADS];
				end[bi+1]=from+lblock[THREADS];
			}
			else
			{
				end[bi+1]=end[bi];
				end[bi]=from+lblock[THREADS];
				beg[bi+1]=to-rblock[THREADS];
			}
			// Increment the stack pointer
			bi++;
		}
		
		__syncthreads();

		unsigned int x = from+lblock[tx+1]-1;
		unsigned int y = to-rblock[tx+1];
		
		coal = from&0xf;

		if(tx+from-coal<to)
		{
			unsigned int d = data[tx+from-coal];

			if(!(tx<coal))
			{
				if(d<pivot)
					data2[x--] = d;
				else
				if(d>pivot)
					data2[y++] = d;
			}
		}

		// Go through the data once again
		// writing it to its correct position
		for(unsigned int i=from+tx+THREADS-coal;i<to;i+=THREADS)
		{	
			unsigned int d = data[i];
			
			if(d<pivot)
				data2[x--] = d;
			else
			if(d>pivot)
				data2[y++] = d;
			
		}

		__syncthreads();

		// As a final step, write the pivot value between the right and left
		// subsequence. Write it to the final destination since this pivot
		// is always correctly sorted
		for(unsigned int i=from+lblock[THREADS]+tx;i<to-rblock[THREADS];i+=THREADS)
		{
			adata[i]=pivot;
		}

		__syncthreads();

	}
#ifdef HASATOMICS
	if(threadIdx.x==0)
		bx = atomicInc(&ohtotal,50000);
	__syncthreads();
#else
	break;
#endif

	}

	__syncthreads();
}
