#include "dtts_merging.incl"

void cgrad (int N, int nnz, float* vals, int* colind, int* rowptr, float* X, float* B, int* niter, float* epsilon){

 sicl_gscsrcg_seq( N, vals, colind, rowptr, X, B,P_NONE,niter,epsilon);
}

void patch_merging(Image* dest_f, Image* patch_f, int dx, int dy, int nlevel, float drange)
{
  patch_merging_cpu(dest_f, patch_f, dx, dy, nlevel, drange);
}
