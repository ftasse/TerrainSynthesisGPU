/**
* GPU Quicksort Testbench
* -----------------------
* Copyright 2007-2008 Daniel Cederman and Philippas Tsigas
*
* This work is licensed under the Creative Commons
* Attribution-Noncommercial-No Derivative Works 3.0
* Unported License. To view a copy of this license,
* visit http://creativecommons.org/licenses/by-nc-nd/3.0/
* or send a letter to Creative Commons, 171 Second Street,
* Suite 300, San Francisco, California, 94105, USA.
*
**/

#define _CRT_SECURE_NO_WARNINGS

#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "gpuqsort.h"
#include "defs.h"
#include "dists.h"

#include <stdlib.h>
#include <string.h>

#ifdef HASSQLITE3
#include "sqlite3.h"
#endif

typedef unsigned int element;


void testDists(char* uniid);
void testModels(char* uniid);
void testPhases(char* uniid);

/**
* The program entrance point
*
* Give no argument to show information
*/
int main(int argc, char* argv[])
{
	if(argc==3&&!strcmp(argv[1],"-d"))
	{
		printf("Testing 'dists'\n");
		testDists(argv[2]);
	}
	else
	if(argc==3&&!strcmp(argv[1],"-m"))
	{
		printf("Testing 'models'\n");
		testModels(argv[2]);
	}
	else
	if(argc==3&&!strcmp(argv[1],"-p"))
	{
		printf("Testing 'phases'\n");
		testPhases(argv[2]);
	}
	else
	{
		printf("\nUsage:\n\ttestbench -[m|d|p] uniqueid\n");
		exit(1);
	}
}


/**
* Checks that data2 contains a sorted list of the elements in data
* @param data   An unsorted list of elements
* @param data2  The possibly correctly sorted version of data
* @param size   The size of the two lists
* @returns      True if data2 is the correctly sorted version of data
*/
bool validate(element* data, element* data2, unsigned int size)
{
	// Sort data using a trusted method
	std::sort(data,data+size);

	// Compare each element to find any differences
	for(unsigned int i=0;i<size;i++)
		if(data[i]!=data2[i])
		{
			// data2 was not correctly sorted
			printf("Error at %d (%i != %i)!\n",i,data[i],data2[i]);
			return false;
		}

		// data2 was correctly sorted
		return true;
}

/**
* Stores the result in a sqlite3 database called "results.sqlite3"
* Creates the file if it does not exist
*/
int saveResults(int testsize, const char* distribution, float time, int threads, int maxblocks, int sbsize, char* uniid, unsigned int phase, unsigned int test)
{
#ifdef HASSQLITE3

	int rc;
	char* zErrMsg;
	char command[2048];
	static sqlite3 *db = 0;

	// Open database if not already open
	if(db==0)
	{
		rc = sqlite3_open("results.sqlite3", &db);
		if(rc)
		{
			fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
			sqlite3_close(db);
			exit(1);
		}

		// Create table to fill in results
		char* createphasetablecmd = "CREATE TABLE results (testsize int, distribution varchar, time float, threads int, maxblocks int, sbsize int, uniid varchar, phase int, test int);";
		sqlite3_exec(db, createphasetablecmd, 0, 0, &zErrMsg);
	}

	// Insert data
	sprintf(command,"insert into results values (%d,\"%s\",%f,%d,%d,%d,\"%s\",%d,%d);\n",testsize, distribution, time, threads, maxblocks, sbsize, uniid, phase, test);
	rc = sqlite3_exec(db, command, 0, 0, &zErrMsg);

	// Exit if we failed to insert data
	if( rc!=SQLITE_OK )
	{
		fprintf(stderr, "SQL error: %s\n", zErrMsg);
		sqlite3_free(zErrMsg);
		sqlite3_close(db);
		exit(1);
	}
#else
	printf("Testsize: %d Distribution: %s Time: %f Threads: %d Maxblocks: %d SBSize: %d uniid: %s phase: %d test: %d\n",testsize, distribution, time, threads, maxblocks, sbsize, uniid, phase, test);
#endif

	return 0;
}

/**
* Tries all distributions ITERATIONS times
*/
void testDists(char* uniid)
{
	const unsigned int MEASURES = 5;
	const unsigned int DISTRIBUTIONS = 6;
	const unsigned int STARTSIZE = 2<<19;

	// Allocate memory for the sequences to be sorted
	unsigned int maxsize = STARTSIZE<<(MEASURES-1);
	element* data = new element[maxsize];
	element* data2 = new element[maxsize];

	double timerValue;
	unsigned int run = 0;

	// Go through all distributions
	for(int d=0;d<DISTRIBUTIONS;d++)
	{
		unsigned int testsize = STARTSIZE;

		// Go through all sizes
		for(int i=0;i<MEASURES;i++,testsize<<=1)
		{
			// Do it several times
			for(int q=0;q<ITERATIONS;q++)
			{
				// Create sequence according to distribution
				dist(data,testsize,d);
				// Store copy of sequence
				memcpy(data2,data,testsize*sizeof(element));

				int threads  =0;
				int maxblocks=0;
				int sbsize   =0;

				// Sort it
				if(gpuqsort(data,testsize,&timerValue,maxblocks,threads,sbsize,0)!=0)
				{
					printf("Error! (%s)\n",getGPUSortErrorStr());
					exit(1);
				}

				// Validate the result
				if(!validate(data2,data,testsize))
					exit(1); 

				saveResults(testsize,getDistString(d),(float)timerValue,threads,maxblocks,sbsize,uniid,0,1);
				printf("%d/%d!\n",run++,MEASURES*DISTRIBUTIONS*ITERATIONS);
			}
		}
	}
}

/**
* Tries all model distributions ITERATIONS times
*/
void testModels(char* uniid)
{
	char* names[] = {"norm-dragon.dat","norm-happy.dat","norm-lucy.dat","norm-manuscript.dat","norm-rgbdragon.dat","norm-statuette.dat"};
	unsigned int run = 0;
	element* data2 = new element[16027872];

	// Try the 6 different models
	for(int i=0;i<6;i++)
	{
		// several times
		for(int q=0;q<ITERATIONS;q++)
		{
			double timerValue;
			int testsize;

			// Read sequence from file
			unsigned int* data = readModel(testsize,names[i]);
			// Store copy for later validation
			memcpy(data2,data,testsize*sizeof(element));

			int threads   =0;
			int maxblocks =0;
			int sbsize    =0;

			// Sort it
			if(gpuqsort(data,testsize,&timerValue,maxblocks,threads,sbsize,0)!=0)
			{
				printf("Error! (%s)\n",getGPUSortErrorStr());
				exit(1);
			}

			// Validate the result
			if(!validate(data2,data,testsize))
				exit(1); 
			saveResults(testsize,names[i],(float)timerValue,threads,maxblocks,sbsize,uniid,0,2);
			printf("%d/%d!\n",run++,6*ITERATIONS);
			free(data);
		}
	}
}

/**
* Tries different combinations of parameters and measures each phase
*/
void testPhases(char* uniid)
{
	unsigned int testsize = 2<<22;
	element* data = new element[testsize];
	element* data2 = new element[testsize];
	double timerValue;

	unsigned int run=0;

	// Uses same distribution for all
	int d = 0;
	dist(data2,testsize,d);

	// Test different sizes
	for(testsize = 2<<18;testsize<=(2<<22)+500;testsize *= 2)  // 5
	{
		// Measure each phase
		for(int phase=0;phase<3;phase++)  // 3
			// Vary the number of threads
			for(int threads=32;threads<=256;threads*=2) // 4
				// Vary the number of blocks
				for(int maxblocks=32;maxblocks<=1024;maxblocks*=2)  // 6
					// Vary when to switch to bitonic
					for(int sbsize=64;sbsize<=2048;sbsize*=2) // 6
						// Do it several times
						for(int q=0;q<ITERATIONS;q++)
						{
							// Store a copy sequence for reuse
							memcpy(data,data2,testsize*4);

							// Sort it
							if(gpuqsort(data,testsize,&timerValue,maxblocks,threads,sbsize,phase)!=0)
							{
								printf("Error! (%s)\n",getGPUSortErrorStr());
								exit(1);
							}

							saveResults(testsize,getDistString(d),(float)timerValue,threads,maxblocks,sbsize,uniid,phase,0);
							printf("%d/%d!\n",run++,2160*ITERATIONS);
						}
	}
}


