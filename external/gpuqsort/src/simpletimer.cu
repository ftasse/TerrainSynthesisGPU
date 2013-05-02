/**
* Simple Timer
* ------------
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


#ifdef _MSC_VER

#include <windows.h>

class SimpleTimer
{
	LARGE_INTEGER starttime;
public:
	void start();
	double end(); 
};

void SimpleTimer::start()
{
	QueryPerformanceCounter(&starttime);
}

double SimpleTimer::end()
{
	LARGE_INTEGER endtime,freq;
	QueryPerformanceCounter(&endtime);
	QueryPerformanceFrequency(&freq);

	return ((double)(endtime.QuadPart-starttime.QuadPart))/((double)(freq.QuadPart/1000.0));
}


#else
#include <sys/time.h>

class SimpleTimer
{
	struct timeval starttime;
public:
	void start();
	double end(); 
};

void SimpleTimer::start()
{
	gettimeofday(&starttime,0);
}

double SimpleTimer::end()
{
	struct timeval endtime;
	gettimeofday(&endtime,0);

	return (endtime.tv_sec - starttime.tv_sec)*1000.0 + (endtime.tv_usec - starttime.tv_usec)/1000.0;
}


#endif
