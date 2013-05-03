#ifndef _INC_NOISE
#define _INC_NOISE
/* file: noise.h
   author: (c) James Gain, 2006
   project: ScapeSketch - sketch-based design of procedural landscapes
   notes: Implementation of Wavlet Noise [Cook & De Rose, "Wavelet Noise", SIGGRAPH2005]
   changes:
*/

#include <cstdlib>
#include <cmath>

#define DEFAULT_LIM 10			// maximum number of noise levels in a hierarchy

#define MATRIX_A 0x9908b0dfUL   // constant vector a
#define UPPER_MASK 0x80000000UL // most significant w-r bits
#define LOWER_MASK 0x7fffffffUL // least significant r bits

using namespace std;

class RandomGen
{
private:

	// Period parameters for random number generator
	static const int N = 624;
	static const int M = 397;

	unsigned long mt[N]; // the array for the state vector
	int mti;

public:

	RandomGen()
	{
		mti = N + 1;
		// initialize random number generator
		unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;
		initByArray(init, length);
	}

	RandomGen(unsigned long init_key[], int key_length)
	{
		mti = N + 1;
		initByArray(init_key, key_length);
	}

	// initGenRand: initializes mt[N] with a seed <s>
	void initGenRand(unsigned long s);

	// initByArray: initialize by an array
	// <init_key> is the array for initializing keys
	// <key_length> is its length
	void initByArray(unsigned long init_key[], int key_length);

	// genRandInt32: generates a random number on [0,0xffffffff]-interval
	unsigned long genRandInt32(void);

	// genRandReal: generates a random number on [0,1)-real-interval
	float genRandReal(void);
};

// useful 1-dimensional wavelet functions
class Wavelet2
{
public:
	// downSample: use uniform B-spline wavelet analysis to reduce a 1D array <from> of <n> coefficients
	//				into an <n>/2 array <to> with <stride> spacing between data points
	//				if <wrap> is true, take cp from opposite side as needed, otherwise repeat the last value on the border
	static void downSample (float *from, float *to, int n, int stride, bool wrap);

	// upSample: use uniform B-spline wavelet synthesis to increase a 1D array <from> of <n>/2 coefficients
	//				into an <n> array <to> with <stride> spacing between data points
	//				if <wrap> is true, take cp from opposite side as needed, otherwise repeat the last value on the border
	static void upSample (float *from, float *to, int n, int stride, bool wrap);

	// upSampleRep: repeat at boundaries rather than wrapping
	static void upSampleRep( float *from, float *to, int n, int stride);
};

// NoiseStats: relevant statistics for a single band of wavelet noise
struct NoiseStats
{
	float min, max, mean, zerodev, meandev, sigma;
	// minimum, maximum, standard deviation from zero
	// the mean is invariably zero for everything except the base level
	// and standard deviation from the mean, and bandwidth (sigma) of a noise signal
};


// NoiseBandStats: collected statistics for a number of noise bands
class NoiseBandStats
{
public:

	NoiseStats * statbands;		// band statistics
	int levels;

	NoiseBandStats()
	{
		statbands = new NoiseStats[DEFAULT_LIM];
		for(int i = 0; i < DEFAULT_LIM; i++)
		{
			statbands[i].min = 0.0f; statbands[i].max = 0.0f; statbands[i].mean = 0.0f;
			statbands[i].zerodev = 0.0f; statbands[i].meandev = 0.0f; statbands[i].sigma = 0.0f;
		}
	}

	NoiseBandStats(int numlevels)
	{
		levels = numlevels;
		statbands = new NoiseStats[levels];
		for(int i = 0; i < levels; i++)
		{
			statbands[i].min = 0.0f; statbands[i].max = 0.0f; statbands[i].mean = 0.0f;
			statbands[i].zerodev = 0.0f; statbands[i].meandev = 0.0f; statbands[i].sigma = 0.0f;
		}
	}

	~NoiseBandStats()
	{
		if(statbands != NULL)
			delete [] statbands;
	}

	// evalStats:	calculate and store the statistics for <band>, a single 1d noise signal
	//				at a certain multiresolution level <lvl> with <size> entries and a bandwidth of <sig>
	//				It makes no difference if this is a linearised 2d array
	void evalStats(int level, int size, float sig, float * band);

	// getStdDev:	return the stadard deviation for the noise at multiresolution level, <level>
	inline float getStdDev(int level){ return statbands[level].zerodev; }

	// load, save: load stats from or save them to a file called <filename>
	bool load(char* filename);
	bool save(char* filename);

	// print: show the current statistics for all levels
	void print();
};


class Noise
{

private:

	float ** tileData;	// pyramid of 2D tiles of sizes 2x2, 4x4, ... , 2^limit x 2^limit
	float * weights;    // band weights assigned to each tile
	int limit;			// highest tile dimension is 2^limit
	RandomGen * rnd;		// random number generator
	Wavelet2 wav;		// to access auxillary functions

	// gaussianNoise: return a noise value in the range [-1;1] with a Gaussian distribution
	float gaussianNoise();

	// generateTile: create a noise tile of dim <n> x <n> by downsampling, upsampling and subtracting Gaussian noise
	void generateTile(int n, float * tile);

	inline int Mod(int x, int n) {int m=x%n; return (m<0) ? m+n : m;}

	// evaluate: evaluate wavelet noise at position <p>[0], <p>[1], where p = [0..<n>,0..<n>] of <tile>
	inline float evalTile(int n, float * tile, float p[2])
	{
		// Non-projected 2D noise
		int i, f[3], c[3], mid[3]; // f, c = filter, noise coeff indices

		float w[2][3], t, result = 0.0f;

		// Evaluate quadratic B-spline basis functions
		for (i = 0; i < 2; i++)
		{
			mid[i]=(int) ceil(p[i]-0.5f);
			t=(float) mid[i]-(p[i]-0.5f);
			w[i][0]=t*t*0.5f;
			w[i][2]=(1-t)*(1-t)*0.5f;
			w[i][1]=1.0f-w[i][0]-w[i][2];
		}

		// Evaluate noise by weighting noise coefficients by basis function values
		for(f[1] = -1; f[1] <= 1; f[1]++)
			for(f[0] = -1; f[0] <= 1; f[0]++)
			{
				float weight = 1.0f;
				for (i = 0; i < 2; i++)
				{
					c[i]=Mod(mid[i]+f[i],n);
					weight*=w[i][f[i]+1];
				}
				result += weight * tile[c[1]*n+c[0]];
			}

		return result;
	}

public:

	// constructor for noise with <lim> number of bands
	Noise(int lim)
	{
		rnd = new RandomGen();
		tileData = NULL;
		generateBands(lim);
	}

	~Noise()
	{
		delete rnd;
		if(tileData != NULL)
		{
			for(int i = 0; i < limit; i++)
				delete tileData[i];
			delete tileData;
		}
	}

	// generateBands: create 2D noise bands from dimension 2 up to 2^<lim>, with each element in [0, 1)
	void generateBands(int lim);

	// weightBands: assign <wghts> array of size <limit> to multiply the values in each band
	// <weights> are not forced into an affine combination, so if this is required it will need
	// to be ensured by the calling function
	void weightBands(float * wghts);

	// evalBands: evaluate 2d noise at position <p>[0], <p>[1] in the range [0, 1] x [0, 1]
	float evalBands(float p[2]);

	// evalWeightedBands:	evaluate 2d noise at position <p>[0], <p>[1] modifying each noise band
	//						according to the <weights> provided for each band. This does not actually
	//						modify the noise itself
	float evalWeightedBands(float * wghts, float p[2]);

	// evalWeightBand: as above but returns the weighted noise of a single band indexed by <k>
	inline float evalWeightBand(float wght, float p[2], int k)
	{
		float v = 0.0f, b;
		float q[2];
		float n = 1.0;
		int i;

		for(i = 0; i < k+1; i++) // n = powf(2.0f, (float) (k+1))
			n *= 2.0f;

		q[0] = p[0] * n;
		q[1] = p[1] * n;

		if(wght != 0.0f)
		{
			b = evalTile(n, tileData[k], q); // band value
			v = (wght * b);
		}

		return v;
	}
};


# endif // _INC_NOISE
