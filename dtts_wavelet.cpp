#include "dtts_wavelet.h"

Image fromVector(float* data, int size)
{

    Image res(size,size);



    for (int i=0; i<size; i++)
        for (int j=0; j<size; j++)
        {
            //cout<<"hi "<<i<<" "<<j<<"\n";
            res.setPixel(i,j, data[i+j*size]);
        }
    res.setMax();

    return res;
}

void Wavelet::transform2d (float *input, float *output, int hsize, int vsize,  int nsteps, int sym_ext=-1)
{
    int j;
    int hLowSize = hsize, hHighSize;
    int vLowSize = vsize, vHighSize;

    // If form of extension unspecified, default to symmetric
    // extensions for symmetrical filters and periodic extensions for
    // asymmetrical filters
    //if (sym_ext == -1)   sym_ext = symmetric;

    float *temp_in = new float [2*npad+max(hsize,vsize)];
    float *temp_out = new float [2*npad+max(hsize,vsize)];

    copy (input, output, hsize*vsize);

    while (nsteps--)
    {
        if ((hLowSize <= 2 || vLowSize <= 2) && sym_ext == 1)
        {
            printf ("Reduce # of transform steps or increase signal size");
            printf ("  or switch to periodic extension");
            printf ("Low pass subband is too small");
        }

        // Do a convolution on the low pass portion of each row
        for (j = 0; j < vLowSize; j++)
        {
            // Copy row j to data array
            copy (output+(j*hsize), temp_in+npad, hLowSize);

            // Convolve with low and high pass filters
            transform_step (temp_in, temp_out, hLowSize, sym_ext);

            // Copy back to image
            copy (temp_out+npad, output+(j*hsize), hLowSize);
        }

        // Now do a convolution on the low pass portion of  each column
        for (j = 0; j < hLowSize; j++)
        {
            // Copy column j to data array
            copy (output+j, hsize, temp_in+npad, vLowSize);

            // Convolve with low and high pass filters
            transform_step (temp_in, temp_out, vLowSize, sym_ext);

            // Copy back to image
            copy (temp_out+npad, output+j, hsize, vLowSize);
        }

        // Now convolve low-pass portion again
        hHighSize = hLowSize/2;
        hLowSize = (hLowSize+1)/2;
        vHighSize = vLowSize/2;
        vLowSize = (vLowSize+1)/2;
    }

    delete [] temp_out;
    delete [] temp_in;

    // Image res = fromVector(output,hsize); res.savePGM("/home/flora/test/filtermidp.pgm",res.maxval); //cout<<"hi\n"; cin.get();

}

/*---------------------------------------------------------------------------*/

void Wavelet::invert2d (float *input, float *output, int hsize, int vsize,
                        int nsteps, int sym_ext=-1)
{
    int i, j;

    // If form of extension unspecified, default to symmetric
    // extensions for symmetrical filters and periodic extensions for
    // asymmetrical filters
    //if (sym_ext == -1) sym_ext = symmetric;

    int *hLowSize = new int [nsteps],
    *hHighSize = new int [nsteps];
    int *vLowSize = new int [nsteps],
    *vHighSize = new int [nsteps];

    hLowSize[0] = (hsize+1)/2;
    hHighSize[0] = hsize/2;
    vLowSize[0] = (vsize+1)/2;
    vHighSize[0] = vsize/2;

    for (i = 1; i < nsteps; i++)
    {
        hLowSize[i] = (hLowSize[i-1]+1)/2;
        hHighSize[i] = hLowSize[i-1]/2;
        vLowSize[i] = (vLowSize[i-1]+1)/2;
        vHighSize[i] = vLowSize[i-1]/2;
    }

    float *temp_in = new float [2*npad+max(hsize,vsize)];
    float *temp_out = new float [2*npad+max(hsize,vsize)];

    copy (input, output, hsize*vsize);

    while (nsteps--)
    {
        // Do a reconstruction for each of the columns
        for (j = 0; j < hLowSize[nsteps]+hHighSize[nsteps]; j++)
        {
            // Copy column j to data array
            copy (output+j, hsize, temp_in+npad,
                  vLowSize[nsteps]+vHighSize[nsteps]);

            // Combine low-pass data (first 1/2^n of signal) with high-pass
            // data (next 1/2^n of signal) to get higher resolution low-pass data
            invert_step (temp_in, temp_out,
                         vLowSize[nsteps]+vHighSize[nsteps], sym_ext);

            // Copy back to image
            copy (temp_out+npad, output+j, hsize,
                  vLowSize[nsteps]+vHighSize[nsteps]);
        }

        // Now do a reconstruction pass for each row
        for (j = 0; j < vLowSize[nsteps]+vHighSize[nsteps]; j++)
        {
            // Copy row j to data array
            copy (output + (j*hsize), temp_in+npad,
                  hLowSize[nsteps]+hHighSize[nsteps]);

            // Combine low-pass data (first 1/2^n of signal) with high-pass
            // data (next 1/2^n of signal) to get higher resolution low-pass data
            invert_step (temp_in, temp_out,
                         hLowSize[nsteps]+hHighSize[nsteps], sym_ext);

            // Copy back to image
            copy (temp_out+npad, output + (j*hsize),
                  hLowSize[nsteps]+hHighSize[nsteps]);
        }
    }

    delete [] hLowSize;
    delete [] hHighSize;
    delete [] vLowSize;
    delete [] vHighSize;

    delete [] temp_in;
    delete [] temp_out;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// input and output are padded with npad values at the beginning and
// at the end

void Wavelet::transform_step (float *input, float *output, int size,
                              int sym_ext)
{
    int i, j;
    // cout<<"Down sampling... "<<size<<" "<<size<<"\n";
    int lowSize = (size+1)/2;
    int left_ext, right_ext;

    if (fsize %2)
    {
        // odd filter length
        left_ext = right_ext = 1;
    }
    else
    {
        left_ext = right_ext = 2;
    }

    //if (sym_ext)    symmetric_extension (input, size, left_ext, right_ext, 1);
    //else
    periodic_extension (input, size);

    //                      coarse  detail
    // xxxxxxxxxxxxxxxx --> HHHHHHHHGGGGGGGG
    for (i = 0; i < lowSize; i++)
    {
        output[npad+i] = 0.0;
        for (j = 0; j < fsize; j++)
        {
            output [npad+i] +=    input[npad + 2*i + first + j] *coeff[j];
        }
    }

    for (i = lowSize; i < size; i++)
    {
        output[npad+i] = 0.0;
        for (j = 0; j < fsize; j++)
        {
            output [npad+i] +=    input[npad + 2*(i-lowSize) + first + j] * coeff[j];
        }
    }

}

/*---------------------------------------------------------------------------*/

void Wavelet::invert_step (float *input, float *output, int size, int sym_ext)
{
    int i, j;
    int left_ext, right_ext, symmetry;
    // amount of low and high pass -- if odd # of values, extra will be
    //   low pass
    int lowSize = (size+1)/2, highSize = size/2;

    symmetry = 1;
    if (fsize % 2 == 0)
    {
        // even length filter -- do (2, X) extension
        left_ext = 2;
    }
    else
    {
        // odd length filter -- do (1, X) extension
        left_ext = 1;
    }

    if (size % 2 == 0)
    {
        // even length signal -- do (X, 2) extension
        right_ext = 2;
    }
    else
    {
        // odd length signal -- do (X, 1) extension
        right_ext = 1;
    }

    int tempsize = 2*npad+lowSize;

    float *temp = new float [2*npad+lowSize];
    for (i = 0; i < lowSize; i++)
    {
        temp[npad+i] = input[npad+i];
    }

    // if (sym_ext)
    //   symmetric_extension (temp, lowSize, left_ext, right_ext, symmetry);
    // else
    periodic_extension (temp, lowSize);

    // coarse  detail
    // HHHHHHHHGGGGGGGG --> xxxxxxxxxxxxxxxx
    for (i = 0; i < 2*npad+size; i++)
        output[i] = 0.0;

    int firstIndex = first;
    int lastIndex = fsize - 1 + firstIndex;

    /*for (i = -lastIndex/2; i <= (size-1-firstIndex)/2; i++)  {
       for (j = 0; j < fsize; j++)  {
    output[npad + 2*i + firstIndex + j] +=  temp[npad+i] * coeff[j];
       }
    }*/

    for (i = -lastIndex/2; i <= (size-1-firstIndex)/2; i++)
    {
        for (j = 0; j < fsize; j++)
        {
            int ind = (npad+i)<0?0:((npad+i)>=tempsize?tempsize-1:(npad+i));   //cout<<ind<<" ";
            int ind2 = npad + 2*i + firstIndex + j;
            ind2 = ind2<0?0:(ind2>=size?size-1:ind2);
            output[ind2] +=  temp[ind] * coeff[j];
        }
    }

    left_ext = 2;

    if (fsize % 2 == 0)
    {
        // even length filters
        right_ext = (size % 2 == 0) ? 2 : 1;
        symmetry = -1;
    }
    else
    {
        // odd length filters
        right_ext = (size % 2 == 0) ? 1 : 2;
        symmetry = 1;
    }

    for (i = 0; i < highSize; i++)
    {
        temp[npad+i] = input[npad+lowSize+i];
    }
    //if (sym_ext)
    //  symmetric_extension (temp, highSize, left_ext, right_ext,
//			  symmetry);
    // else
    periodic_extension (temp, highSize);


    firstIndex = first;
    lastIndex =  fsize - 1 +firstIndex;

    /*for (i = -lastIndex/2; i <= (size-1-firstIndex)/2; i++)  {
       for (j = 0; j < fsize; j++)  {
    output[npad + 2*i + firstIndex + j] +=
      temp[npad+i] * coeff[j];
       }
    }*/
    for (i = -lastIndex/2; i <= (size-1-firstIndex)/2; i++)
    {
        for (j = 0; j < fsize; j++)
        {
            int ind = (npad+i)<0?0:((npad+i)>=tempsize?tempsize-1:(npad+i));   //cout<<ind<<" ";
            int ind2 = npad + 2*i + firstIndex + j;
            ind2 = ind2<0?0:(ind2>=size?size-1:ind2);
            output[ind2] +=  temp[ind] * coeff[j];
        }
    }

    delete [] temp;
}

void Wavelet::periodic_extension (float *output, int size)
{
    int first = npad, last = npad + size-1;

    // extend left periodically
    while (first > 0)
    {
        first--;
        output[first] = output[first+size];
    }

    // extend right periodically
    while (last < 2*npad+size-1)
    {
        last++;
        output[last] = output[last-size];
    }
}

Image Wavelet::down_sample(Image& img,int nsteps)
{
    int hs = img.getHeight();
    int ws = img.getWidth();

    float *input = new float [ws*hs];
    float *output = new float [ws*hs];

    for (int i=0; i<ws; i++)
        for (int j=0; j<hs; j++)
        {
            input[i+j*ws] = img.getPixel(i,j);
            output[i+j*ws] = 0;
        }

    transform2d(input,output,ws,hs,nsteps);

    //Image res = fromVector(output,ws); res.savePGM("/home/flora/test/filtermid.pgm",res.maxval); //cout<<"hi\n"; cin.get();

    Image dest(ws,hs);
    dest.maxval = img.maxval;

    float mini=1000000000.0;
    float maxi=0.;
    for (int k=0; k<ws*hs; k++)
    {
        if (output[k]<mini)  mini=output[k];
        if (output[k]>maxi)  maxi=output[k];
    }

    for (int i=0; i<ws; i++)
        for (int j=0; j<hs; j++)
        {
            dest.setPixel(i,j, (output[i+j*ws]-mini)*img.maxval/maxi);
        }

    delete [] input;
    delete [] output;

    return dest;

}

Image Wavelet::up_sample(Image& img, int nsteps)
{
    int hs = img.getHeight();
    int ws = img.getWidth();

    float *input = new float [ws*hs];
    float *output = new float [ws*hs];

    for (int i=0; i<ws; i++)
        for (int j=0; j<hs; j++)
        {
            input[i+j*ws] = img.getPixel(i,j);
            output[i+j*ws] = 0;
        }

    invert2d(input,output,ws,hs,nsteps);

    //Image res = fromVector(output,ws); res.savePGM("/home/flora/test/filtermid.pgm",res.maxval); //cout<<"hi\n"; cin.get();

    Image dest(ws,hs);
    dest.maxval = img.maxval;

    float mini=1000000000.0;
    float maxi=0.;
    for (int k=0; k<ws*hs; k++)
    {
        if (output[k]<mini)  mini=output[k];
        if (output[k]>maxi)  maxi=output[k];
    }

    for (int i=0; i<ws; i++)
        for (int j=0; j<hs; j++)
        {
            dest.setPixel(i,j, (output[i+j*ws]-mini)*img.maxval/maxi);
        }

    delete [] input;
    delete [] output;

    return dest;

}

Image Wavelet::noise(Image& img, int nsteps)
{
    Image down = down_sample(img,nsteps);
    Image up = up_sample(img,nsteps);
    return diff(up,img);
}


Wavelet::Wavelet()
{
    fsize=0;
    npad=0;
    first = 0;
    coeff=NULL;
}

Wavelet::~Wavelet()
{
    fsize=0;
    npad=0;
    first=0;
    if (coeff!=NULL)    delete [] coeff;
}

Wavelet::Wavelet(const float fcoeff[], int size, int f)
{
    npad=0;
    fsize = size;
    first = f;
    coeff = new float [fsize];
    for (int i=0; i<fsize; i++)
        coeff[i] = fcoeff[i];

}
