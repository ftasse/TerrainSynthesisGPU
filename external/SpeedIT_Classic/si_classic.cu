/*
 * si_classic.cu - part of the SpeedIT Classic toolkit
 * Copyright (C) 2009 - 2010 Vratis
 * email: support@vratis.com
 *
 * SpeedIT Classic toolkit is a free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SpeedIT Classic library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "cuda.h"
#include "cublas.h"
#include <iostream>
#include <limits>

#include "si_classic.h"

using namespace std ;


const int BLOCK_SIZE = 256;

////////////////////////////////////////////////////////////////////////////////
//
//                          Internal functions
//
////////////////////////////////////////////////////////////////////////////////

//------------------------------------------------------------------------------
//
//                              Utilities
//
//------------------------------------------------------------------------------

//
//  Function translates cublas error codes to simple string message.
//
/*static const char*
cublasErrStr(int err_code)
{
  switch(err_code)
  {
    case CUBLAS_STATUS_SUCCESS         :
      return "CUBLAS_STATUS_SUCCESS"          ;
    case CUBLAS_STATUS_NOT_INITIALIZED :
      return "CUBLAS_STATUS_NOT_INITIALIZED"  ;
    case CUBLAS_STATUS_ALLOC_FAILED    :
      return "CUBLAS_STATUS_ALLOC_FAILED"     ;
    case CUBLAS_STATUS_INVALID_VALUE   :
      return "CUBLAS_STATUS_INVALID_VALUE"    ;
    case CUBLAS_STATUS_ARCH_MISMATCH   :
      return "CUBLAS_STATUS_ARCH_MISMATCH"    ;
    case CUBLAS_STATUS_MAPPING_ERROR   :
      return "CUBLAS_STATUS_MAPPING_ERROR"    ;
    case CUBLAS_STATUS_EXECUTION_FAILED:
      return "CUBLAS_STATUS_EXECUTION_FAILED" ;
    case CUBLAS_STATUS_INTERNAL_ERROR  :
      return "CUBLAS_STATUS_INTERNAL_ERROR"   ;
    default                            :
      return "Unknown CUBLAS ERROR"           ;
  } ;
} ;
*/
//
//  Function mainly used to call CUDA routines and check, if everything
//  gone well.
//
/*static inline void
cudaCall(
                cudaError   err,
          const char      * msg
        )
{
  if( cudaSuccess != err)
  {
    cerr << msg << " : "              ;
    cerr << cudaGetErrorString( err ) ;
    cerr << "\n"                      ;
    exit(-1)                          ;
  } ;
} ;

//
//  Function mainly used to call cublas routines and check, if everything
//  gone well.
//
static inline void
cublasCall(
                  cublasStatus   err,
            const char         * msg
          )
{
  if (CUBLAS_STATUS_SUCCESS != err)
  {
    cerr << msg << " : "        ;
    cerr << cublasErrStr( err ) ;
    cerr << "\n"                ;
    exit(-1)                    ;
  } ;
} ;
*/

//------------------------------------------------------------------------------
//
//                            CUDA kernels
//
//------------------------------------------------------------------------------

//
//  GPU kernel to multiply sparse matrix in CSR format by dense vector
//
static __global__ void
KERNEL_crs_multiply(        int     n_rows     ,  // numer of matrix rows
                      const float * vals       ,  // matrix nonzeros
                      const int   * col_idx    ,  // column indices fo nonzeros
                      const int   * row_offset ,  // positions of rows begin, start from 0
                      const float * X          ,  // vector, by which matrix is multiplied
                            float * R             // result
                    )
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (idx >= n_rows)
    return;

  float sum = 0.0 ;

  for (int k = row_offset[ idx ]; k < row_offset[idx + 1]; k++)
  {
    sum +=  vals[ k ]*X[ col_idx[ k ] ];
  }

  R[ idx ] = sum ;
} ;


static __global__ void
KERNEL_replace_with_residual (        int    n_rows , // vector sizes
                               const float * B      , // pointers to vectors
                                     float * R        //
                              )
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x ;

  if (idx >= n_rows)
    return;

  R[idx] = B[idx] - R[idx] ;
} ;

//
//  GPU kernel for simple vector operation
//
//      S = R + alpha * V
//
//  where R, S and V are dense vectors. Kernel is used in CG solver.
//
static __global__ void
KERNEL_smul_vadd(
                        int     n_rows ,
                  const float * R      ,
                  const float * V      ,
                        float   alpha  ,
                        float * S
                )
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x ;

  if (idx >= n_rows) return ;

  S[ idx ] = R[ idx ] + alpha * V[ idx ] ;
} ;


//
//  GPU kernel for simple vector operation
//
//      X = X + alpha * P
//
//  where X and P are dense vectors. Kernel is used in CG solver.
//
static __global__ void
KERNEL_add_to_x (
                        int     n_rows ,
                        float   alpha  ,
                  const float * P      ,
                        float * X
                )
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (idx >= n_rows)  // TODO: tu byl error, poprawic w trunku
    return;

  X[ idx ] += alpha * P[ idx ];
} ;


//------------------------------------------------------------------------------
//
//                          Computational routines
//
//------------------------------------------------------------------------------

//
//  Wrapper for cublasSnrm2() function
//
static inline float
norm(
      const float * v,
      const int     size
    )
{
  float result = cublasSnrm2(size, v, 1)      ;

  return result;
} ;

//
//  Wrapper for cublasSdot() function
//
static inline float
dot_product (
                    int n_rows ,
              const float* v   ,
              const float* w
            )
{
  float result =  cublasSdot(n_rows, v, 1, w, 1)     ;

  return result ;
} ;

//
//  Function used in CG solver
//
static void
calc_residual(      int     n_rows ,
              const float * vals   ,
              const int   * c_idx  ,
              const int   * r_idx  ,
              const float * B      ,
              const float * X      ,
                    float * R
              )
{
  const int size = n_rows ;
  int numBlocks = size / BLOCK_SIZE + (size % BLOCK_SIZE == 0 ? 0 : 1) ;
  dim3 dimGrid(numBlocks) ;
  dim3 dimBlock(BLOCK_SIZE) ;

  //cout<<"Kernel size: "<<n_rows<<" "<<dimBlock.x*dimGrid.x<<endl;
  sicl_gscsrmv(n_rows, vals, c_idx, r_idx, X, R);  // R = A*X;

  // R = B - R
  KERNEL_replace_with_residual <<<dimGrid, dimBlock>>> (n_rows, B, R) ;

  //cudaThreadSynchronize() ;
  //cudaCall (cudaGetLastError(), "KERNEL_replace_with_residual FAILED") ;
} ;


//
//  S = R + alpha * V
//
static void
smul_vadd(
                  int     n_rows ,
            const float * R      ,
            const float * V      ,
                  float   alpha  ,
                  float * S
         )
{
  const int size = n_rows;
  int numBlocks = size / BLOCK_SIZE + (size % BLOCK_SIZE == 0 ? 0 : 1) ;
  dim3 dimGrid(numBlocks) ;
  dim3 dimBlock(BLOCK_SIZE) ;

  KERNEL_smul_vadd<<<dimGrid,dimBlock>>>( n_rows, R, V, alpha, S) ;

  //cudaThreadSynchronize() ;
  //cudaCall (cudaGetLastError(), "KERNEL_smul_vadd FAILED") ;
} ;

//
//  X = X + alpha * P
//
static void
add_to_x(
                int     n_rows ,
                float   alpha  ,
          const float * P      ,
                float * X
        )
{
  const int size = n_rows ;
  int numBlocks = size / BLOCK_SIZE + (size % BLOCK_SIZE == 0 ? 0 : 1);
  dim3 dimGrid(numBlocks);
  dim3 dimBlock(BLOCK_SIZE);

  KERNEL_add_to_x <<<dimGrid,dimBlock>>> (n_rows, alpha, P, X);

  //cudaThreadSynchronize();
  //cudaCall (cudaGetLastError(), "KERNEL_add_to_x FAILED") ;
} ;

////////////////////////////////////////////////////////////////////////////////
//
//                            Exported functions
//
////////////////////////////////////////////////////////////////////////////////


//
//  Function computes sparse matrix by dense vector multiplication
//
//            A * X = R
//
//  where A - matrix, X and R - vectors. Computation is done with
//  GPU.
//
int
sicl_gscsrmv(
                int     n_rows ,  // numer of matrix rows
          const float * vals   ,  // matrix nonzeros
          const int   * c_idx  ,  // column indices fo nonzeros
          const int   * r_idx  ,  // positions of rows begin, start from 0
          const float * X      ,  // vector, by which matrix is multiplied
                float * R         // result
            )
{
  const int size = n_rows;
  int numBlocks = size / BLOCK_SIZE + (size % BLOCK_SIZE == 0 ? 0 : 1);
  dim3 dimGrid(numBlocks);
  dim3 dimBlock(BLOCK_SIZE);

  //cout<<"Kernel launch: "<<n_rows<<" "<<dimBlock.x*dimGrid.x<<endl;
  //cout<<"Launch kernel: "<<size<<" "<<dimGrid.x*dimBlock.x<<endl;

  KERNEL_crs_multiply<<<dimGrid,dimBlock>>> (
                                              n_rows ,
                                              vals   ,
                                              c_idx  ,
                                              r_idx  ,
                                              X      ,
                                              R
                                            ) ;
  //cudaThreadSynchronize();
  //cudaCall (cudaGetLastError(), "KERNEL_crs_multiply FAILED") ;

  return 0 ;
} ;

//
//  Function solves with preconditioned Conjugate Gradient (CG) method system of
//  linear algebraic equations
//
//      A * X = B
//
//  where A is a sparse matrix and X and B are dense vectors.
//
int
sicl_gscsrcg(     int            n_rows  ,
            const float        * vals    ,
            const int          * c_idx   ,
            const int          * r_idx   ,
                  float        * X       ,
            const float        * B       ,
                  PRECOND_TYPE   precond , // Only for compatibility
                  int          * n_iter  ,
                  float        * epsilon
             )
{
  int result = -1 ;
  const float almost_zero = numeric_limits<float>::min();

  // Notation, algorithm: see Barlett et al, p. 13
  float rho_1;      // \rho_{i-1}
  float rho_2 = 0;  // \rho_{i-2}
  float alpha = 0;  // \alpha{i}
  float beta;       // \beta_{i-1}

  float norm_b = norm(B, n_rows) ;

  if (norm_b < almost_zero) // if B == 0
  {
    norm_b = 1.0 ;
  }

  float *R, *P, *Z, *Q ;

  // residuals
  cudaMalloc ((void**)(&R), n_rows*sizeof(float)) ;
  // search directions
  cudaMalloc ((void**)(&P), n_rows*sizeof(float));
  // solution of A * Z = R
  cudaMalloc ((void**)(&Z), n_rows*sizeof(float));
  cudaMalloc ((void**)(&Q), n_rows*sizeof(float));

  float* pZ = R ;
  //cudaCall(cudaMalloc ((void**)(&pZ), n_rows*sizeof(float)), "cudaMalloc failed for pZ") ;
  //cudaCall(cudaMemcpy(pZ, R, n_rows*sizeof(float), cudaMemcpyDeviceToDevice),          "cudaMemcpy pZ = R FAILED" ) ;

  // R = B - AX;  (Line 1 at Barlett's algorithm)
  calc_residual(n_rows, vals, c_idx, r_idx, B, X, R) ;

  float residuum = norm(R, n_rows) / norm_b ;

  if (residuum < *epsilon)  // if the trial solution satisfies the equation...
  {
    *n_iter = 0 ;
    result  = 0 ;
  }
  else
  {
    for (int niter = 1 ; niter <= *n_iter ; niter++)
    {
      rho_1 = dot_product(n_rows, R, pZ) ;

      if (1 == niter) {
        // p^1 = z^0;   Barlett: line 6
        cudaMemcpy(P, pZ, n_rows*sizeof(float), cudaMemcpyDeviceToDevice);
      } else {
        beta = (rho_1/rho_2) ;
        smul_vadd(n_rows, pZ, P, beta, P) ;             // P = Z + beta * P
      } ;

      sicl_gscsrmv(n_rows, vals, c_idx, r_idx, P, Q) ; // Q = A*P
      alpha = rho_1 / dot_product(n_rows, P, Q)  ;

      add_to_x(n_rows,  alpha, P, X) ;         // X^i = X^{i-1} + alpha_i * p^i
      add_to_x(n_rows, -alpha, Q, R) ;
      rho_2 = rho_1 ;

      residuum = norm(R, n_rows) / norm_b ;

      if (residuum < *epsilon) // iteration succeeded
      {
        *n_iter = niter ;
        result  = 0 ;
        break ;
      } ;
    } ;
  } ;

  *epsilon = residuum ;

  cudaFree (Q);
  cudaFree (R);
  cudaFree (P);
  cudaFree (Z);


  return result ;
} ;
