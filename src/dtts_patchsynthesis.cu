#include "dtts_patchsynthesis.incl"

#define WIDTH 800
#define HEIGHT 600


void terrain_synthesis(char* exemplar, char* sketchf, char* output, char ridges, int bsize)
{

   	int num_devices, device;
	cudaGetDeviceCount(&num_devices);
	if (num_devices > 1) {
	      int max_multiprocessors = 0, max_device = 0;
	      for (device = 0; device < num_devices; device++) {
		      cudaDeviceProp properties;
		      cudaGetDeviceProperties(&properties, device);
		      if (max_multiprocessors < properties.multiProcessorCount) {
		              max_multiprocessors = properties.multiProcessorCount;
		              max_device = device;
		      }
	      }
	      cudaSetDevice(max_device);
	}
	
	terrain_synthesis_run(exemplar, sketchf, output, ridges, bsize);
}
