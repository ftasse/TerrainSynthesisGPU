/* file: noise.cpp
   author: (c) James Gain, 2006
   project: ScapeSketch - sketch-based design of procedural landscapes
   notes: Implementation of Wavlet Noise [Cook & De Rose, "Wavelet Noise", SIGGRAPH2005]
   changes:
*/

#include "noise.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>

/*
   Modified from A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved

   DISTRIBUTION NOTICE REMOVED - will need to be added back in if
   distributed
*/

///////////////////
//// RANDOMGEN ////
///////////////////

void RandomGen::initGenRand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++)
	{
        mt[mti] =
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
        // See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier.
        // In the previous versions, MSBs of the seed affect
        // only MSBs of the array mt[]

        mt[mti] &= 0xffffffffUL;
        // for >32 bit machines
    }
}

void RandomGen::initByArray(unsigned long init_key[], int key_length)
{
    int i, j, k;
    initGenRand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--)
	{
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL)) + init_key[j] + j; // non linear
        mt[i] &= 0xffffffffUL; // for WORDSIZE > 32 machines
        i++; j++;
        if (i>=N)
		{
			mt[0] = mt[N-1];
			i=1;
		}
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--)
	{
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL)) - i; // non linear
        mt[i] &= 0xffffffffUL; // for WORDSIZE > 32 machined
        i++;
        if (i>=N)
		{
			mt[0] = mt[N-1]; i=1;
		}
    }

    mt[0] = 0x80000000UL; // MSB is 1; assuring non-zero initial array
}

unsigned long RandomGen::genRandInt32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    // mag01[x] = x * MATRIX_A  for x=0,1

    if (mti >= N)
	{ // generate N words at one time
        int kk;

        if (mti == N+1)   // if init_genrand() has not been called,
            initGenRand(5489UL); // a default initial seed is used

        for (kk=0;kk<N-M;kk++)
		{
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++)
		{
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }

    y = mt[mti++];

    // Tempering
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

// genRandReal: generates a random number on [0,1)-real-interval
float RandomGen::genRandReal(void)
{
    return genRandInt32()*(1.0/4294967296.0);
    // divided by 2^32
}


///////////////
//// NOISE ////
///////////////

float Noise::gaussianNoise()
{
	float x1, x2;
	float ret, r2;

	do {
		x1 = 2.0f * (float) rnd->genRandReal() - 1.0f;	// [-1, 1)
		x2 = 2.0f * (float) rnd->genRandReal() - 1.0f;

		r2 = x1*x1 + x2*x2;
	} while ((r2 == 0) || (r2 > 1.0));

	ret = x1 * sqrtf((-2.0f * logf(r2))/r2);
	ret *= 0.25f;		// Possibility of ( N(0, 1) < 4.0 ) = 100%

	if (ret < -1.0) ret = -1.0; // Account for loss of precision
	if (ret >  1.0) ret = 1.0;

	return ret;
}

// Note: this code is designed for brevity, not efficiency; many operations can be hoisted,
// precomputed, or vectorized. Some of the straightforward details, such as tile meshing,
// decorrelating bands and fading out the last band, are omitted in the interest of space.

int Mod(int x, int n) {int m=x%n; return (m<0) ? m+n : m;}
bool inRange(int x, int n) {return (x >= 0 && x < n);}

#define ARAD 16

void Wavelet2::downSample (float *from, float *to, int n, int stride, bool wrap)
{
	int i, k;
	float *a, aCoeffs[2*ARAD] = {
				0.000334f,-0.001528f, 0.000410f, 0.003545f,-0.000938f,-0.008233f, 0.002172f, 0.019120f,
				-0.005040f,-0.044412f, 0.011655f, 0.103311f,-0.025936f,-0.243780f, 0.033979f, 0.655340f,
				0.655340f, 0.033979f,-0.243780f,-0.025936f, 0.103311f, 0.011655f,-0.044412f,-0.005040f,
				0.019120f, 0.002172f,-0.008233f,-0.000938f, 0.003546f, 0.000410f,-0.001528f, 0.000334f};

	a = &aCoeffs[ARAD];
	for (i=0; i<n/2; i++)
	{
		to[i*stride] = 0.0f;
		for (k=2*i-ARAD; k<2*i+ARAD; k++)
		{
			if(wrap)
				to[i*stride] += a[k-2*i] * from[Mod(k,n)*stride];
			else
				if(inRange(k, n))
					to[i*stride] += a[k-2*i] * from[k*stride];
		}
	}
}

void Wavelet2::upSample( float *from, float *to, int n, int stride, bool wrap)
{
	int i, k;
	float *p, pCoeffs[4] = { 0.25f, 0.75f, 0.75f, 0.25f };
	p = &pCoeffs[2];
	for (i=0; i<n; i++)
	{
		to[i*stride] = 0.0f;
		for (k=i/2; k<=i/2+1; k++)
		{
			if(wrap)
				to[i*stride] += p[i-2*k] * from[Mod(k,n/2)*stride];
			else
				if(inRange(k, n/2))
					to[i*stride] += p[i-2*k] * from[k*stride];
		}
	}
}

void Wavelet2::upSampleRep( float *from, float *to, int n, int stride)
{
	int i, k;
	float *p, pCoeffs[4] = { 0.25f, 0.75f, 0.75f, 0.25f };
	p = &pCoeffs[2];
	float strfrom=0.;

	for (i=0; i<n; i++)
	{
		to[i*stride] = 0.0f;
		for (k=i/2; k<=i/2+1; k++)
		{
			if(inRange(k, n/2))
			{
				to[i*stride] += p[i-2*k] * from[k*stride];
				strfrom = from[k*stride];
			}
			else
				to[i*stride] += p[i-2*k] * strfrom;

		}
	}
}

void Noise::generateTile(int n, float * tile)
{
	if (n%2) n++; // tile size must be even

	int ix, iy, i, sz=n*n;
	float * temp1 = new float[sz];
	float * temp2 = new float[sz];

	// Step 1. Fill the tile with random numbers in the range -1 to 1.
	for (i=0; i<n*n; i++)
		tile[i] = gaussianNoise();

	// Steps 2 and 3. Downsample and upsample the tile
	for (iy=0; iy<n; iy++) // each y row
	{
		i = iy*n;
		wav.downSample( &tile[i], &temp1[i], n, 1, true);
		wav.upSample( &temp1[i], &temp2[i], n, 1, true);

	}

	for (ix=0; ix<n; ix++) // each x row
	{
		i = ix;
		wav.downSample( &temp2[i], &temp1[i], n, n, true);
		wav.upSample( &temp1[i], &temp2[i], n, n, true);
	}

	// Step 4. Subtract out the coarse-scale contribution
	for (i=0; i<n*n; i++)
		tile[i]-=temp2[i];

	/*
	int j;
	bool outofrange = false;

	for(j = 0; j < n*n; j++)
		if(noise[j] > 1.0f || noise[j] < -1.0f)
			outofrange = true;

	if(outofrange)
		wxLogMessage(wxT("Noise is out of [-1, 1] range"));
	*/

	// map noise onto [0, 1) range
	//for(i = 0; i < n*n; i++)
	//	tile[i] = (tile[i] + 0.2f); // / 2.0f;

	// Avoid even/odd variance difference by adding odd-offset version of noise to itself
	//	but this messes with the range so avoid here
//	int offset=n/2;
//	if (offset%2==0)
//		offset++;
//	for (i=0,ix=0; ix<n; ix++)
//		for (iy=0; iy<n; iy++)
//			temp1[i++] = noise[ Mod(ix+offset,n) + Mod(iy+offset,n)*n ];

//	for (i=0; i<n*n; i++)
//		noise[i]+=temp1[i];

	delete temp1;
	delete temp2;
}

void Noise::generateBands(int lim)
{
	int n = 2;
	limit = lim;

	tileData = new float *[limit];
	weights = new float[limit];

	// generate tiles as increasing powers of 2
	for(int i = 0; i < limit; i++)
	{
		tileData[i] = new float[n*n];
		generateTile(n, tileData[i]);
		n *= 2;
	}
}

void Noise::weightBands(float * wghts)
{
	int i;
	float sum = 0.0f;

	// assign to internal weighting array
	for(i = 0; i < limit; i++)
	{
		weights[i] = wghts[i];
		sum += weights[i];
	}

	/*
	// normalize if necessary
	if(sum < 1.0f - pluszero || sum > 1.0f + pluszero)
		for(i = 0; i < limit; i++)
			weights[i] /= sum;
	*/
}

float Noise::evalBands(float p[2])
{
	float v = 0.0;
	float q[2];
	int n = 2;

	for(int i = 0; i < limit; i++)
	{
		q[0] = p[0] * (float) n;
		q[1] = p[1] * (float) n;
		if(weights[i] > 0.0f) // only eval if necessary
			v += weights[i] * evalTile(n, tileData[i], q);
		n *= 2;
	}

	return v;
}

float Noise::evalWeightedBands(float * wghts, float p[2])
{
	float v = 0.0f, b;
	float q[2];
	int n = 2;

	for(int i = 0; i < limit; i++)
	{
		q[0] = p[0] * (float) n;
		q[1] = p[1] * (float) n;

		if(wghts[i] != 0.0f)
		{
			b = evalTile(n, tileData[i], q); // band value
			v += (wghts[i] * b);
		}
		n *= 2;
	}

	return v;
}

////////////////////////
//// NoiseBandStats ////
////////////////////////

void NoiseBandStats::evalStats(int level, int size, float sig, float * band)
{
	float cmax, cmin, csum, cdev;
	int n, j;

	//n = (int) powf(2.0f, (float) (level+1));
	n = size;
	cmax = -100.0f; cmin = 100.0f; csum = 0.0f;
	statbands[level].sigma = sig;

	for(j = 0; j < n; j++)
	{
		if(band[j] > cmax)
			cmax = band[j];
		if(band[j] < cmin)
			cmin = band[j];
		csum += band[j];
	}
	statbands[level].min = cmin;
	statbands[level].max = cmax;
	statbands[level].mean = csum / (float) n;

	// calculate standard deviation
	// assuming that the mean is zero
	csum = 0.0f;
	for(j = 0; j < n; j++)
	{
		cdev = band[j]; // - statbands[level].mean;
		cdev *= cdev;
		csum += cdev;
	}
	statbands[level].zerodev = sqrtf(csum / (float) (n-1));

	// now calculate the standard deviation from the mean
	csum = 0.0f;
	for(j = 0; j < n; j++)
	{
		cdev = band[j] - statbands[level].mean;
		cdev *= cdev;
		csum += cdev;
	}
	statbands[level].meandev = sqrtf(csum / (float) (n-1));
}


bool NoiseBandStats::load(char* filename)
{

	ifstream infile;

	infile.open((char *) filename, ios_base::in);
	if(infile.is_open() && !infile.eof())
	{
		infile >> levels;
		delete statbands;
		statbands = new NoiseStats[levels];

		for(int i = 0; i < levels; i++)
			infile >> statbands[i].min;
		for(int i = 0; i < levels; i++)
			infile >> statbands[i].max;
		for(int i = 0; i < levels; i++)
			infile >> statbands[i].mean;
		for(int i = 0; i < levels; i++)
			infile >> statbands[i].zerodev;
		for(int i = 0; i < levels; i++)
			infile >> statbands[i].meandev;
		for(int i = 0; i < levels; i++)
			infile >> statbands[i].sigma;

		infile.close();
		return true;
	}
	else
	{
		infile.close();
		return false;
	}

	return false;
}

bool NoiseBandStats::save(char* filename)
{
	ofstream outfile;
	int i;

	outfile.open((char *) filename, ios_base::out);
	if(outfile.is_open())
	{
		outfile << levels << endl;
		for(i = 0; i < levels; i++) // min
			outfile << statbands[i].min << " ";
		outfile << endl;
		for(i = 0; i < levels; i++) // max
			outfile << statbands[i].max << " ";
		outfile << endl;
		for(i = 0; i < levels; i++) // mean
			outfile << statbands[i].mean << " ";
		outfile << endl;
		for(i = 0; i < levels; i++) // zerodev
			outfile << statbands[i].zerodev << " ";
		outfile << endl;
		for(i = 0; i < levels; i++) // meandev
			outfile << statbands[i].meandev << " ";
		outfile << endl;
		for(i = 0; i < levels; i++) // sigma
			outfile << statbands[i].sigma << " ";
		outfile << endl;
		outfile.close();
		return true;
	}
	else
		return false;
}

// print: show the current statistics for all levels
void NoiseBandStats::print()
{

	// mac

	cout << "num levels =" << levels<<endl;
	for(int i = 0; i < levels; i++)
	{
		cout << "Statistics for Band " << i;
		cout<<endl;
		cout << "  min = "<< statbands[i].min << ", max =" << statbands[i].max;
		cout<<endl;
		cout << " mean = " << statbands[i].mean << ", meandev = " << statbands[i].meandev;
		cout<<endl;
		cout << "  zerodev = " << statbands[i].zerodev << ", sigma = " <<  statbands[i].sigma;
		cout<<endl;
	}

	// windows
/*
	wxLogMessage("num levels = %d", levels);
	for(int i = 0; i < levels; i++)
	{
		wxLogMessage("Statistics for Band %d", i);


		wxLogMessage(wxT("  min = %f, max = %f", statbands[i].min, statbands[i].max));
		wxLogMessage(wxT("  mean = %f, meandev = %f", statbands[i].mean, statbands[i].meandev));
		wxLogMessage(wxT("  zerodev = %f, sigma = %f", statbands[i].zerodev, statbands[i].sigma));
	}
	*/
}
