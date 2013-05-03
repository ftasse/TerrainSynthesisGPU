#if defined(__CUDACC__)
#include "cublas.h"
#include <cusparse.h>
#include <cublas.h>
#endif

#include "dtts_merging.incl"

#define TILE_WIDTH 16

#if defined(__CUDACC__)
__global__
void sinterpolate_gpu(point_t* points,seaminfo* seaminf, int nsize, int bsize, float drange){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    int j = blockIdx.y*blockDim.y + threadIdx.y;

    if (i<bsize && j<bsize && points[i+j*bsize].x>-1e3 ){

    		    node_t p(i,j);
                    float wtot=0, wsrcx=0, wsrcy=0;

                    for (unsigned int k=0; k<nsize; k++)
                    {
			float d  = sqrtf((p.x-seaminf[k].pt.x)*(p.x-seaminf[k].pt.x) + (p.y-seaminf[k].pt.y)*(p.y-seaminf[k].pt.y));

                        float w = powf((drange-d)/(drange*d),4);	//powf(((d*d)/(drange*drange))-1.,2);
                        wtot = wtot+w;
                        if (d <=  drange)
                        {
                            wsrcx += w*(seaminf[k].diffx);
                            wsrcy += w*(seaminf[k].diffy);
                        }

                    }


                    points[i+j*bsize].x = (wsrcx/wtot);
                    points[i+j*bsize].y = (wsrcy/wtot);
    }

}
#endif

#if defined(__CUDACC__)
void cgrad (int N, int nnz, float* vals, int* colind, int* rowptr, float* X, float* B, int* niter, float* epsilon){

 cublasInit();
  //for (int k=0; k<N; k++){
    //   X[k] = 0.0;
 	//cout<<b[k]<<" ";
	//}
//cout<<endl;

 float* vals_dev;        cublasAlloc(nnz, sizeof(float), (void**) &vals_dev);
 int* colind_dev;        cublasAlloc(nnz, sizeof(int), (void**) &colind_dev);
 int * rowptr_dev;       cublasAlloc(N+1, sizeof(int), (void**) &rowptr_dev);
 float * X_dev;       cublasAlloc(N, sizeof(float), (void**) &X_dev);
 float * B_dev;       cublasAlloc(N, sizeof(float), (void**) &B_dev);
 //int* niter_dev;         cublasAlloc(1, sizeof(int), (void**) &niter_dev);
 //float* epsilon_dev;     cublasAlloc(1, sizeof(float), (void**) &epsilon_dev);

  cublasSetVector (nnz, sizeof(float),vals, 1, vals_dev, 1);
  cublasSetVector (nnz, sizeof(int),colind, 1, colind_dev, 1);
  cublasSetVector (N+1, sizeof(int),rowptr, 1, rowptr_dev, 1);
  cublasSetVector (N, sizeof(float),X, 1, X_dev, 1);
  cublasSetVector (N, sizeof(float),B, 1, B_dev, 1);
  //*niter = 0;


/*
cudaDeviceProp deviceProp;
    int devID = 0;
    if (devID < 0) {
       printf("exiting...\n");
       exit(0);
    }
    cudaGetDeviceProperties(&deviceProp, devID) ;
 printf("> GPU device has %d Multi-Processors, SM %d.%d compute capabilities\n\n",
		deviceProp.multiProcessorCount, deviceProp.major, deviceProp.minor);

    int version = (deviceProp.major * 0x10 + deviceProp.minor);
    if(version < 0x11)
    {
        printf("Requires a minimum CUDA compute 1.1 capability\n");
        printf("PASSED");
        cudaThreadExit();
    }*/

  //sicl_gscsrcg_seq( N, vals, colind, rowptr, X, B,P_NONE,niter,epsilon);
  sicl_gscsrcg( N, vals_dev, colind_dev, rowptr_dev, X_dev, B_dev,P_NONE,niter,epsilon);
  //bicgstab_kernel( N, vals_dev, colind_dev, rowptr_dev, X_dev, B_dev,P_NONE,niter,epsilon);
  //sicl_gscsrmv( N, vals_dev, colind_dev, rowptr_dev, X_dev, B_dev);

  /*int max_iter  =10000;

    cusparseHandle_t handle = 0;
    cusparseStatus_t status;
    status = cusparseCreate(&handle);
    if (status != CUSPARSE_STATUS_SUCCESS) {
        fprintf( stderr, "!!!! CUSPARSE initialization error\n" );
        return ;
    }

    cusparseMatDescr_t descr = 0;
    status = cusparseCreateMatDescr(&descr);
    if (status != CUSPARSE_STATUS_SUCCESS) {
        fprintf( stderr, "!!!! CUSPARSE cusparseCreateMatDescr error\n" );
        return ;
    }
    cusparseSetMatType(descr,CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(descr,CUSPARSE_INDEX_BASE_ZERO);

    float a, b, r0, r1;
    float *d_Ax;
    float *d_p;
    cudaMalloc((void**)&d_p, N*sizeof(float));
    cudaMalloc((void**)&d_Ax, N*sizeof(float));

    cusparseScsrmv(handle,CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, 1.0, descr, vals_dev, rowptr_dev, colind_dev, X_dev, 0.0, d_Ax);
    cublasSaxpy(N, -1.0, d_Ax, 1, B_dev, 1);
    r1 = cublasSdot(N, B_dev, 1, B_dev, 1);

    int k = 1;
    const float tol = 1e-5;

    while (r1 > tol*tol && k <= max_iter) {
        if (k > 1) {
            b = r1 / r0;
            cublasSscal(N, b, d_p, 1);
            cublasSaxpy(N, 1.0, B_dev, 1, d_p, 1);
        } else {
            cublasScopy(N, B_dev, 1, d_p, 1);
        }

        cusparseScsrmv(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, 1.0, descr, vals_dev, rowptr_dev, colind_dev,d_p, 0.0, d_Ax);
        a = r1 / cublasSdot(N, d_p, 1, d_Ax, 1);
        cublasSaxpy(N, a, d_p, 1, X_dev, 1);
        cublasSaxpy(N, -a, d_Ax, 1, B_dev, 1);

        r0 = r1;
        r1 = cublasSdot(N, B_dev, 1, B_dev, 1);
        cudaThreadSynchronize();
        //shrLog("iteration = %3d, residual = %e\n", k, sqrtf(r1));
        k++;
    }

 cudaFree(d_p);
    cudaFree(d_Ax);
*/ cublasGetVector (N, sizeof(float),X_dev, 1, X, 1);

  cublasFree(vals_dev);
  cublasFree(colind_dev);
  cublasFree(rowptr_dev);
  cublasFree(X_dev);
  cublasFree(B_dev);
  //cublasFree(niter_dev);
  //cublasFree(epsilon_dev);

  cublasShutdown();

  //cout<<"Niter: "<<*niter<<" Epsilon: "<<*epsilon<<endl;
  //for (int k=0; k<N; k++)	cout<<X[k]<<" ";
  //cout<<endl;
}

#if defined(__CUDACC__)
void patch_merging(Image* dest_f, Image* patch_f, int dx, int dy, int nlevel, float drange)
{
    int w = (*patch_f).width(), h = patch_f->height();

    Image* mask = graphCut(dest_f,patch_f,dx,dy);    //Tested and severe=false performs the best
    //mask->savePGM("/tmp/mask.pgm");

    Gradient pdest_g(Image(w,h),Image(w,h));
    Gradient patch_g = get_gradient(*patch_f);
    vector<node_t> seam;
    vector<float>   xdiff, ydiff;

    vector<seaminfo> seaminf;

     for (int i=0; i<w; i++)
        for (int j=0; j<h; j++){
            if ( (*mask)(i,j) >= vsSINK){
                float h = dest_f->getPixelXY(i+dx, j+dy);
                float h1 = dest_f->getPixelXY(i+dx-1, j+dy);
                float h2 = dest_f->getPixelXY(i+dx, j+dy-1);
                if (h1>BG) (pdest_g.first)(i,j)= h - h1;
                if (h2>BG) (pdest_g.second)(i,j)= h - h2;
            }
        }

    for (int j=0; j<patch_f->height(); j++)
		   for (int i=0; i<patch_f->width(); i++){
            if ( (*mask)(i,j) == vsSINK ){
                //seam.push_back(node_t(i,j));
                //xdiff.push_back(pdest_g.first(i,j)-patch_g.first(i,j));
                //ydiff.push_back(pdest_g.second(i,j)-patch_g.second(i,j));

                seam.push_back(node_t(i,j));
                seaminfo tmp(node_t(i,j), pdest_g.first(i,j)-patch_g.first(i,j), pdest_g.second(i,j)-patch_g.second(i,j) );
		seaminf.push_back(tmp);
            }

            if ( (*mask)(i,j) >= vsSINK ){
                patch_g.first(i,j) = pdest_g.first(i,j);
                patch_g.second(i,j) = pdest_g.second(i,j);
                (*patch_f)(i,j) = dest_f->getPixelXY(i+dx, j+dy);
            }

        }

    //patch_g.first.savePGM("/tmp/patch_gx_baf.pgm",20);    patch_g.second.savePGM("/tmp/patch_gy_baf.pgm",20);
    //pdest_g.first.savePGM("/tmp/patch_gx_bcf.pgm",20);    pdest_g.second.savePGM("/tmp/patch_gy_bcf.pgm",20);

    if (seam.size()>0)
    {

	{
		int bsize = mask->width();
		//patch_g.first.savePGM("/tmp/patch_gx_bef.pgm",20);    patch_g.second.savePGM("/tmp/patch_gy_bef.pgm",20);
		// Shepard Interpolation
		point_t*  points = new point_t[bsize*bsize];
		for (int j=0; j<patch_f->height(); j++)
		   for (int i=0; i<patch_f->width(); i++){
				 point_t tmp(-1e5,-1e5);
				 if ( (*mask)(i,j) < vsSINK )
				 	tmp = point_t(0,0);
				 points[i+j*bsize] = tmp;
			}

		//sinterpolate_cpu(points, &seaminf[0], seaminf.size(), bsize, drange);

		point_t* points_dev;	cudaMalloc( (void**) &points_dev, sizeof(point_t)*bsize*bsize );
		seaminfo* seaminf_dev;	cudaMalloc( (void**) &seaminf_dev, sizeof(seaminfo)*seaminf.size() );

		cudaMemcpy(points_dev, points,  sizeof(point_t)*bsize*bsize, cudaMemcpyHostToDevice );
		cudaMemcpy(seaminf_dev, &seaminf[0],  sizeof(seaminfo)*seaminf.size(), cudaMemcpyHostToDevice );

		dim3 dimGrid( (bsize/TILE_WIDTH)+1, (bsize/TILE_WIDTH)+1);
	   	dim3 dimBlock(TILE_WIDTH,TILE_WIDTH);
		sinterpolate_gpu<<<dimGrid,dimBlock>>>(points_dev, seaminf_dev, seaminf.size(), bsize, drange);
		cudaMemcpy(points, points_dev,  sizeof(point_t)*bsize*bsize, cudaMemcpyDeviceToHost );

		cudaFree(points_dev);
		cudaFree(seaminf_dev);

		for (int j=0; j<patch_f->height(); j++)
		   for (int i=0; i<patch_f->width(); i++)
		    	if ( (*mask)(i,j) < vsSINK ){
		    	patch_g.first(i,j) = patch_g.first(i,j) + points[i+j*bsize].x; //(wsrcx/wtot);
		        patch_g.second(i,j) = patch_g.second(i,j) + points[i+j*bsize].y; //(wsrcy/wtot);
		    }

		delete [] points;
	}

        //patch_g.first.savePGM("/tmp/patch_gx_bef.pgm",20);    patch_g.second.savePGM("/tmp/patch_gy_bef.pgm",20);
        //sinterpolate_g(patch_g,*mask,seam,xdiff,ydiff,drange);
        //patch_g.first.savePGM("/tmp/patch_gx.pgm",20);    patch_g.second.savePGM("/tmp/patch_gy.pgm",20);

        Image div= get_divergent(patch_g);
        //div.savePGM("/tmp/patch_div.pgm",15);

        int *pos = new int [w*h];
        uint N = 0;

        for (int y=0; y<h; y++)
            for (int x=0; x<w; x++)
            {
                pos[x+y*w] = N;
                N++;
            }
        //cin.get();
        //patch_f->savePGM("/tmp/res_cand1.pgm",dest_f->maxval);
        poissonsolve(dest_f,patch_f,div,pos,N,dx,dy);
        //patch_f->savePGM("/tmp/res_cand2.pgm",dest_f->maxval);
        //cin.get();
        delete [] pos;

    }


    for (int i=0; i<patch_f->width(); i++)
        for (int j=0; j<patch_f->height(); j++)
        {
            if (dest_f->inBounds(i+dx,j+dy))
            {
                (*dest_f)( i+dx,j+dy) = (*patch_f)(i,j) ;
            }
        }

    delete mask;
}
#endif
