#ifndef DTTS_MERGING_H
#define DTTS_MERGING_H

#include "external/maxflow/graph.h"

#include "dtts_image.h"
using namespace Dtts;

#define MAX_CAPACITY MAX_FLOAT
typedef Graph<float,float,float> GraphType;
typedef pair<node_t,bool> cut_node;

typedef pair<Image,Image> Gradient;

Gradient get_gradient(Image& src);

Gradient get_gradient_r(Image& src);

Image get_divergent(Gradient grad);

Image get_divergent(Image& src);

float w_eqn(float x);

int near2seam( node_t p, vector<node_t> seam);

void wire_deform_shepard(Image pdest, Image& psrc, Image& mask, vector<node_t> seam, Image& target, float doff);

float graphCut_cost(Image* dest, Image* patch, int dx, int dy, bool gradient = true, bool severe = false);

Image* graphCut(Image* dest, Image* patch, int dx, int dy, bool gradient = true, bool severe = false);

void paste_cut (Image* dest, Image* patch, int dx, int dy);

void poissonsolve(Image* dest, Gradient& grad, float boundary=0);

void poisson(Image* dest, Image* patch, int dx, int dy, int nlevel=3, float drange=25);

void poisson_blend(Image* dest, Image* patch, int dx, int dy, int nlevel=3, float drange=25);

void shepard(Image* dest, Image* patch, int dx, int dy, int nlevel=3, float drange=25);

void wiredef(Image* dest, Image* patch, int dx, int dy, int nlevel=3, float doff=25);

void patch_merging(Image* dest_f, Image* patch_f, int dx, int dy, int nlevel=3, float drange=25);
void patch_merging_cpu(Image* dest_f, Image* patch_f, int dx, int dy, int nlevel=3, float drange=25);

void feathering(Image* dest, Image* patch, int dx, int dy, int nlevel=3, float drange=25);

#endif
