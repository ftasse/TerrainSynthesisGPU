#ifndef DTTS_IMAGE_H
#define DTTS_IMAGE_H

#pragma once

#include "headers.h"

/* Seed the Random number generator */
void initrand();

/* generates a psuedo-random float between 0.0 and 0.999... */
float randfloat();

/* generates a psuedo-random float between 0.0 and max */
float randfloat(float max);

/* generates a psuedo-random float between min and max */
float randfloat(float min, float max);


const float laplace_filter[][3] =  {	{0,-1,0}, {-1, 4,-1},  {0, -1, 0}  };
const float sobel_filterX[][3] = {	{-1, -2, -1}, {0, 0, 0},  {1, 2, 1}  } ;
const float sobel_filterY[][3] = {	{-1, 0, 1}, {-2, 0, 2},  {-1, 0, 1}  } ;

extern float ndistance(node_t p, node_t q);
extern float mindist(node_t pnode, node_list nodes);

bool node_in(const node_t pnode, const node_list plist);
vector<node_t> get_rotate_pts(vector<node_t> nodes, node_t ctr, int theta);
void rotate_pts(vector<node_t>& nodes, node_t ctr, int theta);
void mirrorX_pts(vector<node_t>& nodes, int bsize);
void mirrorY_pts(vector<node_t>& nodes, int bsize);


/*---------------------------------------------------------------------------*/
// helpful inline functions
/*---------------------------------------------------------------------------*/
#define rmin(x,y) (((x)<(y))?(x):(y))
#define rmax(x,y) (((x)>(y))?(x):(y))

/*---------------------------------------------------------------------------*/
inline float mod (float x, float N) {
   float xmodN = x - N*((int)(x/N));
   if (xmodN < 0) xmodN += N;
   return xmodN;
}

/*---------------------------------------------------------------------------*/
inline float square (float x) { return (x*x); }
/*---------------------------------------------------------------------------*/
inline int  isquare (int x) { return (x*x); }
/*---------------------------------------------------------------------------*/
inline int  sign (float x)   { return (x > 0 ? 1 : x < 0 ? -1 : 0); }
/*---------------------------------------------------------------------------*/
inline int log2 (int x) {
   int count = 0;

   while (x > 1)  {
      x >>= 1;
      count++;
   }
   return count;
}

// Dtts name space
namespace Dtts {

    /**
     * Image class. Simple image in texture synthesis :-)
     *
     * @author Flora Tasse
     * @version 1.0.0
     */
    class Image{

    public:

        /**
         * Construct an Image instance
         */
        Image();

        /**
         * Construct an Image instance with a given size
         *@param pwidth new width of the image
         *@param pheight new height of the image
         */
        Image(int pwidth, int pheight);

        /**
         * Destuction of an Image instance
         */
        ~Image();

        /**
         * Load a PGM image
         *
         *@param fname path of the PGM image
         */
        void loadPGM(const char* fname);

        /**
         * Save to a PGM image
         *
         *@param fname path of the PGM image
         */
        void savePGM(const char* fname, size_t threshold=MAX_VAL);

        void convertTerragen(const char* fname);

        /**
         * Set the width of the image
         *
         *@param pwidth new width of the image
         */
        void setWidth(const int pwidth);

        /**
         * Set the height of the image
         *
         *@param pheight new height of the image
         */
        void setHeight(const int pheight);


        /**
         * Set the size of the image
         *
         *@param pwidth new width of the image
         *@param pheight new height of the image
         */
        void setSize(const int pwidth, const int pheight);

         /**
         * Set the value of a pixel
         *
         *@param i column of the pixel
         *@param j row of the pixel
         *@param pval new value of the pixel at position (i,j)
         */
        void setPixel(const int i, const int j, const float pval);

        /**
         * Set the value of a pixel with Neumann boundary conditions
         *
         *@param i column of the pixel
         *@param j row of the pixel
         *@param pval new value of the pixel at position (i,j)
         */
        void setPixelXY(const int i, const int j, const float pval);

         /**
         * Get the width of the image
         *
         *@return width of the image in pixels
         */
         int getWidth();

         /**
         * Get the height of the image
         *
         *@return heigh of the image in pixels
         */
        int getHeight();

        /**
         * Get the value of a pixel
         *
         *@param i column of the pixel
         *@param j row of the pixel
         *
         *@return pixel value at position (i,j)
         */
        float getPixel(const int i, const int j);

        /**
         * Get the value of a pixel with Neumann boundary conditions
         *
         *@param i column of the pixel
         *@param j row of the pixel
         *
         *@return pixel value at position (i,j)
         */
        float getPixelXY(const int i, const int j);
        float linearXY(const float i, const float j);
        float cubicXY(const float i, const float j);

        /**
         * Get the value of a pixel
         *
         *@param pnode position of the pixel
         *
         *@return pixel value at position pnode
         */
        float getPixel(const node_t pnode);

        /**
         * Get the value of a pixel
         *
         *@param pnode position of the pixel
         *
         *@return pixel value at position pnode
         */
        float getPixelXY(const node_t pnode);

        /**
         * Blur the image using a gaussian kernel
         */
        void gaussianBlur();

        /**
         * Remove detail from the image using wavelet compression
         */
        void compress();

        void reverse();

        Image(Image& img,int rate);

        vector<node_t> on_line(node_t p, node_t q);

        void drawline(node_t p, node_t q, float color);

        void setMax();
        float getMax();
        float getMin();

        int width(){   return mwidth; }
        int height(){   return mheight; }
        float atXY(const node_t pnode) {   return getPixelXY(pnode); }
        float atXY(const int i, const int j)  {    return getPixelXY(i,j); }
        float& operator()(const int i, const int j=0) { return mpixels[i+j*mwidth]; }

        Image get_crop(int x0, int y0, int x1, int y1);
        Image get_crop(int x0, int y0, int x1, int y1,float theta);
        Image get_crop(int x0, int y0, int x1, int y1,float theta, int mir);

        void resize(int nwidth, int nheight);

        float* getPixels(){    return mpixels; }

        Image filter();
	    void filter_me();
        Image filter(int nlevels);
        Image get_rotate(int theta);
        void rotate(int theta);

        Image convolute(const float kernel[][3], int ksize);

        vector<Image> get_pyramid(int nlevels);

        Image(const Image& oldImage);
        void operator=(const Image&);
        bool inBounds(int i, int j);
        Image (Image& img, int x, int y, int width,int height);

        float maxval;
        float minval;

        void mirrorX();
        void mirrorY();
        void mirror(int m);

        void normalize(){
            float vmax = maxval;
            for (int j=0; j<height(); j++)
                for (int i=0; i<width(); i++){
                    float nval = mpixels[i+j*mwidth]*255./vmax;
                    mpixels[i+j*mwidth] = nval;
                }

        }
        void unnormalize(){
            float vmax = maxval;
            for (int j=0; j<height(); j++)
                for (int i=0; i<width(); i++){
                    float nval = mpixels[i+j*mwidth]*vmax/255.;
                    mpixels[i+j*mwidth] = nval;
                }

        }

    protected:

        // Width of the image
        int mwidth;

        // Height of the image
        int mheight;

        // Pixel values of the image
        float* mpixels;

    };

     /**
      * Terrain class. Simple terrain derived from Image :-)
      *
      * @author Flora Tasse
      * @version 1.0.0
      */
     class Terrain : public Image {

     public:

        Terrain();

        /**
         * Load a Terragen file
         *
         *@param fname path of the .ter file
         */
        bool loadTerragen(const char* fname);

        //Terrain(const Terrain& oldImage);
        void operator=(const Terrain&);

        /**
         * Save to a Terragen file
         *
         *@param fname path of the .ter file
         */
        bool saveTerragen(const char* fname);

        /**
         * Set the elevation at a given point
         *
         *@param i horizontal position
         *@param j vertical position
         *@param eval new elevation at position (i,j)
         */
        void setElevation(const int i, const int j, const int eval);

        /**
         * Get the elevation at a given point
         *
         *@param i horizontal position
         *@param j vertical position
         *
         *@return elevation at position (i,j) in terrain units
         */
        int getElevation(const int i, const int j);

        /**
         *  Render the terrain using points
         */
        //void render(float terscale, float hscale, GLenum rendermode);

        // Scale of the terrain in meters per terrain units
        float mscale;

        // Base height of the terrain
        int mbaseheight;

       // Scale of the height values of the terrain
        int mheightscale;

        //GLuint terrainDL;

        Image getImage();

	void operator=(const Image&);

	Terrain(int pwidth, int pheight);

    private:



        /* ****** Rendering ******** */



    };



}

Dtts::Image sobelEdge(Dtts::Image);
Dtts::Image diff(Dtts::Image& img1,Dtts::Image& img2);
Dtts::Image add(Dtts::Image& img1,Dtts::Image& img2);

float   compare_pnsr(Dtts::Image& img1, Dtts::Image&img2 );  // Compares one image to another and returns the PSNR
float   compare_mse(Dtts::Image& img1, Dtts::Image&img2 );  // Compares one image to another and returns the mean squared error
float   compare_overlap(Dtts::Image& img1, Dtts::Image&img2, int x, int y );

vector<point_t> points_on_line(point_t p, point_t q, int nsteps);

#if defined (__CUDACC__)
__host__ __device__
#endif
inline void swap(int & a, int & b)
{
    // Alternative swap doesn't use a temporary register:
    a ^= b;
    b ^= a;
    a ^= b;
}

void normalizeImage(Dtts::Image &image );

#endif
