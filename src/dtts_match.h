#include "dtts_image.h"
using namespace Dtts;

#include "dtts_ppa.h"
#include "dtts_wavelet.h"
#include "dtts_2dmorph.h"
#include "dtts_merging.h"
#include "external/noise/noise.h"

#define KNUM 5
#define DROT   45
#define DMIR  2

#define use_cut 0
#define use_ssd 1
#define use_bend 0
#define use_profile 1

#define use_angle 1
#define use_noisestat 1


void paste_patch(Image& dest, Image& patch, int drange, int dx, int dy);

void match_noFeature(Terrain& dest, Image& src, Image& target, node_list dem_nodes, vector<Image>& tar_pyr, vector<Image>&  src_pyr, int bsize, int osize);

void match_Feature(Terrain& dest, Tree& usr_features, Tree& dem_features, vector<Image>& tar_pyr, vector<Image>&  src_pyr, int bsize);

//void match_noFeature_cpu(Terrain& dest, Image& src, Image& target, node_list dem_nodes, vector<Image>& tar_pyr, vector<Image>&  src_pyr, int bsize, int osize);

//void match_Feature_cpu(Terrain& dest, Tree& usr_features, Tree& dem_features, vector<Image>& tar_pyr, vector<Image>&  src_pyr, int bsize);

void patch_synthesis(Terrain& dest, Terrain& src,int bsize=100, int osize=25);

void sort_host(cost_t* values, int num);
