/**
* GPU Quicksort Library
* ---------------------
* Copyright 2007-2008 Daniel Cederman and Philippas Tsigas
*
* This work is licensed under the Creative Commons
* Attribution-Noncommercial-No Derivative Works 3.0
* Unported License. To view a copy of this license,
* visit http://creativecommons.org/licenses/by-nc-nd/3.0/
* or send a letter to Creative Commons, 171 Second Street,
* Suite 300, San Francisco, California, 94105, USA.
**/

#define BUILDING_DLL

#include "stdio.h"
#include "gpuqsort.h"

#include "simpletimer.cu"

#include <algorithm>

// Keep tracks of the data blocks in phase one
template <typename element>
struct BlockSize
{
	unsigned int beg;
	unsigned int end;
	unsigned int orgbeg;
	unsigned int orgend;
	element		 rmaxpiv;
	element		 lmaxpiv;
	element		 rminpiv;
	element		 lminpiv;

	bool		 altered;
	bool		 flip;
	element		 pivot;
};

// Holds parameters to the kernel in phase one
template <typename element>
struct Params
{
	unsigned int from;
	unsigned int end;
	element pivot;
	unsigned int ptr;
	bool last;
};

// Used to perform a cumulative sum between blocks.
// Unnecessary for cards with atomic operations.
// Will be removed when these becomes more common
template <typename element>
struct Length
{
	element maxpiv[MAXBLOCKS];
	element minpiv[MAXBLOCKS];

	unsigned int left[MAXBLOCKS];
	unsigned int right[MAXBLOCKS];
};

// Since we have divided up the kernel in to three
// we need to remember the result of the cumulative sum
// Unnecessary for cards with atomic operations.
// Will be removed when these becomes more common
struct Hist
{
	unsigned int left[(MAXTHREADS)*MAXBLOCKS];
	unsigned int right[(MAXTHREADS)*MAXBLOCKS];
};

struct LQSortParams
{
	unsigned int beg;
	unsigned int end;
	bool flip;
	unsigned int sbsize;
};

#include "gpuqsort_kernels.cu"

#undef THREADS
#define THREADS threads

/**
* The main sort function
* @param data		Data to be sorted
* @param size		The length of the data
* @param timerValue Contains the time it took to sort the data [Optional]
* @returns 0 if successful. For non-zero values, use getErrorStr() for more information about why it failed.
*/
template <typename element>
int GPUQSort<element>::sort(element* data, unsigned int size, double* timerValue, unsigned int blockscount, unsigned int threads, unsigned int sbsize, unsigned int phase)
{
	if(!init)
		return 1;

	if(!threads||!blockscount||!sbsize)
	{
		threads   = 1<<(int)round(log(size * TK + TM)/log(2.0));
		blockscount = 1<<(int)round(log(size * MK + MM)/log(2.0));
		sbsize    = 1<<(int)round(log(size * SK + SM)/log(2.0));
	}

#ifdef HASATOMICS
		unsigned int* doh;
		unsigned int oh;

		cudaGetSymbolAddress((void**)&doh,"ohtotal");
		oh=0;
		cudaMemcpy(doh,&oh,4,cudaMemcpyHostToDevice);
#endif

	if(threads>MAXTHREADS)
		return 1; 

	if(blockscount>MAXBLOCKS)
		return 1;

	SimpleTimer st;

	// Copy the data to the graphics card and create an auxiallary array
	ddata2 = 0; ddata = 0;
	if(!errCheck(cudaMalloc((void**)&ddata2,(size)*sizeof(element))))
		return 1;
	if(!errCheck(cudaMalloc((void**)&ddata,(size)*sizeof(element))))
		return 1;
	if(!errCheck(cudaMemcpy(ddata, data, size*sizeof(element), cudaMemcpyHostToDevice) ))
		return 1;


	if(timerValue!=0)
	{
		// Start measuring time
		cudaThreadSynchronize();
		
		st.start();
	}

	// We start with a set containg only the sequence to be sorted
	// This will grow as we partition the data
	workset[0].beg = 0;
	workset[0].end = size;
	workset[0].orgbeg = 0;
	workset[0].orgend = size;
	workset[0].altered = false;
	workset[0].flip = false;

	// Get a starting pivot
	workset[0].pivot = (min(min(data[0],data[size/2]),data[size-1]) + max(max(data[0],data[size/2]),data[size-1]))/2;
	unsigned int worksize = 1;

	unsigned int blocks = blockscount/2;
	unsigned totsize = size;
	unsigned int maxlength = (size/blocks)/4;

	unsigned int iterations = 0;
	bool flip = true;

	// Partition the sequences until we have enough
	while(worksize<blocks)
	{
		unsigned int ws = totsize/blocks;
		unsigned int paramsize = 0;

		// Go through the sequences we have and divide them into sections
		// and assign thread blocks according to their size
		for(unsigned int i=0;i<worksize;i++)
		{
			if((workset[i].end-workset[i].beg)<maxlength)
				continue;

			// Larger sequences gets more thread blocks assigned to them
			unsigned int blocksassigned = max((workset[i].end-workset[i].beg)/ws,1);
			for(unsigned int q=0;q<blocksassigned;q++)
			{
				params[paramsize].from = workset[i].beg + ws*q;
				params[paramsize].end = params[paramsize].from + ws;
				params[paramsize].pivot = workset[i].pivot;
				params[paramsize].ptr = i;
				params[paramsize].last = false;
				paramsize++;
				
			}
			params[paramsize-1].last = true;
			params[paramsize-1].end = workset[i].end;

			workset[i].lmaxpiv=0;
			workset[i].lminpiv=0xffffffff;
			workset[i].rmaxpiv=0;
			workset[i].rminpiv=0xffffffff;
		}

		if(paramsize==0)
			break;

		// Copy the block assignment to the GPU
		if(!errCheck(cudaMemcpy(dparams, params, paramsize*sizeof(Params<element>), cudaMemcpyHostToDevice) ))
			return 1;

		// Do the cumulative sum
		if(flip)
			part1<<< paramsize, THREADS, (THREADS+1)*2*4+THREADS*2*4 >>>(ddata,dparams,dhists,dlength);
		else
			part1<<< paramsize, THREADS, (THREADS+1)*2*4+THREADS*2*4 >>>(ddata2,dparams,dhists,dlength);
		if(!errCheck((cudaMemcpy(length, dlength,sizeof(Length<element>) , cudaMemcpyDeviceToHost) )))
			return 1; 

		// Do the block cumulative sum. Done on the CPU since not all cards have support for
		// atomic operations yet. 
		for(unsigned int i=0;i<paramsize;i++)
		{
			unsigned int l = length->left[i];
			unsigned int r = length->right[i];
			
			length->left[i] = workset[params[i].ptr].beg;
			length->right[i] = workset[params[i].ptr].end;
			
			workset[params[i].ptr].beg+=l;
			workset[params[i].ptr].end-=r;
			workset[params[i].ptr].altered = true;
			
			workset[params[i].ptr].rmaxpiv = max(length->maxpiv[i],workset[params[i].ptr].rmaxpiv);
			workset[params[i].ptr].lminpiv = min(length->minpiv[i],workset[params[i].ptr].lminpiv);
			
			workset[params[i].ptr].lmaxpiv = min(workset[params[i].ptr].pivot,workset[params[i].ptr].rmaxpiv); 
			workset[params[i].ptr].rminpiv = max(workset[params[i].ptr].pivot,workset[params[i].ptr].lminpiv); 

			
		}

		// Copy the result of the block cumulative sum to the GPU
		if(!errCheck((cudaMemcpy(dlength, length, sizeof(Length<element>), cudaMemcpyHostToDevice) )))
			return 1;

		// Move the elements to their correct position
		if(flip)
			part2<<< paramsize, THREADS >>>(ddata,ddata2,dparams,dhists,dlength);
		else
			part2<<< paramsize, THREADS >>>(ddata2,ddata,dparams,dhists,dlength);

		// Fill in the pivot value between the left and right blocks
		part3<<< paramsize, THREADS >>>(ddata,dparams,dhists,dlength);

		flip = !flip;

		// Add the sequences resulting from the partitioning
		// to set
		unsigned int oldworksize = worksize;
		totsize = 0;
		for(unsigned int i=0;i<oldworksize;i++)
		{
			if(workset[i].altered)
			{
				if(workset[i].beg-workset[i].orgbeg>=maxlength)
					totsize += workset[i].beg-workset[i].orgbeg;
				if(workset[i].orgend-workset[i].end>=maxlength)
					totsize += workset[i].orgend-workset[i].end;

				workset[worksize].beg = workset[worksize].orgbeg = workset[i].orgbeg;
				workset[worksize].end = workset[worksize].orgend = workset[i].beg;
				workset[worksize].flip=flip;
				workset[worksize].altered = false;
				workset[worksize].pivot = (workset[i].lminpiv/2+workset[i].lmaxpiv/2);

				worksize++;

				workset[i].orgbeg = workset[i].beg = workset[i].end;
				workset[i].end = workset[i].orgend;
				workset[i].flip=flip;
				workset[i].pivot = (workset[i].rminpiv/2+workset[i].rmaxpiv/2);
				workset[i].altered = false;
			}
		}
		iterations++;

	}

	// Due to the poor scheduler on some graphics card
	// we need to sort the order in which the blocks
	// are sorted to avoid poor scheduling decisions
	unsigned int sortblocks[MAXBLOCKS*2];
	for(int i=0;i<worksize;i++)
		sortblocks[i]=((workset[i].end-workset[i].beg)<<(int)round(log((float)(MAXBLOCKS*4.0f))/log(2.0f))) + i;
	std::sort(&sortblocks[0],&sortblocks[worksize]);

	if(worksize!=0)
	{
		// Copy the block assignments to the GPU
		for(int i=0;i<worksize;i++)
		{
		 	unsigned int q = (worksize-1)-(sortblocks[i]&(MAXBLOCKS*4-1));

			lqparams[i].beg =  workset[q].beg;
			lqparams[i].end = workset[q].end;
			lqparams[i].flip = workset[q].flip;
			lqparams[i].sbsize = sbsize;
		}

		if(!errCheck((cudaMemcpy(dlqparams, lqparams, worksize*sizeof(LQSortParams), cudaMemcpyHostToDevice) )))
			return 1;

		// Run the local quicksort, the one that doesn't need inter-block synchronization
		if(phase!=1) 
			lqsort<<< worksize, THREADS, max((THREADS+1)*2*4,sbsize*4) >>>(ddata,ddata2,dlqparams,phase);
	}

	cudaThreadSynchronize();

	if(timerValue!=0)
	{
		// Measure the time taken
		*timerValue = st.end();
	}

	err = cudaThreadSynchronize();
	// Free the data
	if(err!=cudaSuccess)
	{
		cudaFree(ddata);
		cudaFree(ddata2);
		return 1;
	}

	// Copy the result back to the CPU
	if(!errCheck((cudaMemcpy(data, ddata, size*sizeof(element), cudaMemcpyDeviceToHost) )))
		return 1;

	cudaFree(ddata);
	cudaFree(ddata2);

	return 0;
}

template <typename element>
bool GPUQSort<element>::errCheck(int e)
{
	if(e==cudaSuccess)
		return true;

	err = e;
	cudaFree(ddata);
	cudaFree(ddata2);
	return false;
}

template <typename element>
GPUQSort<element>::GPUQSort():init(false),workset(0),params(0),length(0),lqparams(0),dlqparams(0),
							  dhists(0),dlength(0),dparams(0)
{
	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, 0);
	if(!strcmp(deviceProp.name,"GeForce 8800 GTX"))
	{
		TK = 1.17125033316e-005f;
		TM = 52.855721393f;
		MK = 3.7480010661e-005f;
		MM = 476.338308458f;
		SK = 4.68500133262e-005f;
		SM = 211.422885572f;
	}
	else
	if(!strcmp(deviceProp.name,"GeForce 8600 GTS"))
	{
		TK = 0.0f;
		TM = 64.0f;
		MK = 0.0000951623403898f;
		MM = 476.338308458f;
		SK = 0.0000321583081317f;
		SM = 202.666666667f;
	}
	else
	{
		TK = 0;
		TM = 128;
		MK = 0;
		MM = 512;
		SK = 0;
		SM = 512;
	}

	if(cudaMallocHost((void**)&workset,MAXBLOCKS*2*sizeof(BlockSize<element>))!=cudaSuccess) return;
	if(cudaMallocHost((void**)&params,MAXBLOCKS*sizeof(Params<element>))!=cudaSuccess) return;
	if(cudaMallocHost((void**)&length,sizeof(Length<element>))!=cudaSuccess) return;
	if(cudaMallocHost((void**)&lqparams,MAXBLOCKS*sizeof(LQSortParams))!=cudaSuccess) return;
	if(cudaMalloc((void**)&dlqparams,MAXBLOCKS*sizeof(LQSortParams))!=cudaSuccess) return;
	if(cudaMalloc((void**)&dhists,sizeof(Hist))!=cudaSuccess) return;
	if(cudaMalloc((void**)&dlength,sizeof(Length<element>))!=cudaSuccess) return;
	if(cudaMalloc((void**)&dparams,MAXBLOCKS*sizeof(Params<element>))!=cudaSuccess) return;

	init = true;
}

/**
* Returns the latest error message
* @returns the latest error message
*/
template <typename element>
const char* GPUQSort<element>::getErrorStr()
{
	return cudaGetErrorString((cudaError_t)err);
}

template <typename element>
GPUQSort<element>::~GPUQSort()
{
	cudaFreeHost(workset);
	cudaFreeHost(params);
	cudaFreeHost(length);
	cudaFreeHost(lqparams);
	cudaFree(dparams);
	cudaFree(dlqparams);
	cudaFree(dhists);
	cudaFree(dlength);
}

// Exported functions

char* expErrMsg = "No errors";

 GPUQSort<unsigned int>* s=0;

extern "C" 
DLLEXPORT int gpuqsort(unsigned int* data, unsigned int size, double* timerValue, unsigned int blockscount, unsigned int threads, unsigned int sbsize, unsigned int phase)
{
	if(s==0)
		s=new GPUQSort<unsigned int>();


	if(s->sort(data,size,timerValue, blockscount, threads, sbsize, phase)!=0)
	{
		expErrMsg = (char*)s->getErrorStr();
		return 1;
	}
	else
		return 0;
}

// Float support removed due to some problems with CUDA 2.0 and templates
// Will be fixed

//extern "C" DLLEXPORT 
/*int gpuqsortf(float* data, unsigned int size, double* timerValue)
{
	GPUQSort<float> s;
	if(s.sort(data,size,timerValue)!=0)
	{
		expErrMsg = (char*)s.getErrorStr();
		return 1;
	}
	else
		return 0;
}*/

extern "C"
DLLEXPORT const char* getGPUSortErrorStr()
{
	return expErrMsg;
}
