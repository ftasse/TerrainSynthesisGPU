#define Sqrt2 1.414213562

#include "dtts_image.h"
using namespace Dtts;

const float HaarCoeffs [] = { 1.0/Sqrt2, 1.0/Sqrt2 };

// A few Daubechies filters

const float Daub4Coeffs [] = { 0.4829629131445341,  0.8365163037378077,
                                0.2241438680420134, -0.1294095225512603
                              };

const float Daub6Coeffs [] = { 0.3326705529500825,  0.8068915093110924,
                                0.4598775021184914, -0.1350110200102546,
                                -0.0854412738820267,  0.0352262918857095
                              };

const float Daub8Coeffs [] = { 0.2303778133088964,  0.7148465705529154,
                                0.6308807679398587, -0.0279837694168599,
                                -0.1870348117190931,  0.0308413818355607,
                                0.0328830116668852, -0.0105974017850690
                              };


class Wavelet
{

    void transform_step (float *input, float *output, int size, int sym_ext);
    void invert_step (float *input, float *output, int size, int sym_ext);


public:

    int fsize, npad, first;
    float *coeff;

    Wavelet();
    ~Wavelet();
    Wavelet(const float fcoeff[], int size, int f);

    void transform2d (float *input, float *output,  int hsize, int vsize,  int nsteps, int sym_ext);
    void invert2d (float *input, float *output, int hsize, int vsize,  int nsteps, int sym_ext);

    void periodic_extension (float *output, int size);

    Image down_sample(Image& dest, int nsteps);
    Image up_sample(Image& dest, int nsteps);
    Image noise(Image& img, int nsteps);


    //copy length elements from p1 to p2
    void copy (const float *p1, float *p2, const int length)
    {
        int temp = length;
        while(temp--) *p2++ = *p1++;
    }
    void copy (const float *p1, const int stride1, float *p2, const int length)
    {
        int temp = length;
        while(temp--)
        {
            *p2++ = *p1;
            p1 += stride1;
        }
    }
    void copy (const float *p1, float *p2, const int stride2, const int length)
    {
        int temp = length;
        while(temp--)
        {
            *p2 = *p1++;
            p2 += stride2;
        }
    }

};
