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

#include <iostream>
#include <limits>

#include "si_classic.h"
#include <cmath>

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

//------------------------------------------------------------------------------
//
//                            CUDA kernels
//
//------------------------------------------------------------------------------

//
//  GPU kernel to multiply sparse matrix in CSR format by dense vector
//

//
//  CPU kernel to multiply sparse matrix in CSR format by dense vector
//
static  void
KERNEL_crs_multiply_seq(        int     n_rows     ,  // numer of matrix rows
                      const float * vals       ,  // matrix nonzeros
                      const int   * col_idx    ,  // column indices fo nonzeros
                      const int   * row_offset ,  // positions of rows begin, start from 0
                      const float * X          ,  // vector, by which matrix is multiplied
                            float * R             // result
                    )
{
  for (int idx=0; idx<n_rows; idx++){

	  float sum = 0.0 ;

	  for (int k = row_offset[ idx ]; k < row_offset[idx + 1]; k++)
	  {
	    sum +=  vals[ k ]*X[ col_idx[ k ] ];
	  }

	  R[ idx ] = sum ;
  }

} ;

//
//  GPU kernel for simple vector operation
//
//      R = B - R
//
//  where R and B are dense vectors. Kernel is used in CG solver.
//

static void
KERNEL_replace_with_residual_seq (        int    n_rows , // vector sizes
                               const float * B      , // pointers to vectors
                                     float * R        //
                              )
{
   for (int idx=0; idx<n_rows; idx++){

	 R[idx] = B[idx] - R[idx] ;
  }
} ;

static void
KERNEL_smul_vadd_seq(
                        int     n_rows ,
                  const float * R      ,
                  const float * V      ,
                        float   alpha  ,
                        float * S
                )
{
  for (int idx=0; idx<n_rows; idx++){

 	S[ idx ] = R[ idx ] + alpha * V[ idx ] ;

  }
} ;

//
//  GPU kernel for simple vector operation
//
//      X = X + alpha * P
//
//  where X and P are dense vectors. Kernel is used in CG solver.
//

static void
KERNEL_add_to_x_seq (
                        int     n_rows ,
                        float   alpha  ,
                  const float * P      ,
                        float * X
                )
{
  for (int idx=0; idx<n_rows; idx++){

  X[ idx ] += alpha * P[ idx ];

  }
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
norm_seq(
      const float * v,
      const int     size
    )
{
  float result=0;
  for (int i=0; i<size; i++){
  	result+=v[i]*v[i];
  }
  result = sqrtf(result);
  return result;
} ;

//
//  Wrapper for cublasSdot() function
//

static inline float
dot_product_seq (
                    int n_rows ,
              const float* v   ,
              const float* w
            )
{
	  float result=0;
	  for (int i=0; i<n_rows; i++){
	  	result+=v[i]*w[i];
	  }
	  return result;
} ;

//
//  Function used in CG solver
//

static void
calc_residual_seq(      int     n_rows ,
              const float * vals   ,
              const int   * c_idx  ,
              const int   * r_idx  ,
              const float * B      ,
              const float * X      ,
                    float * R
              )
{
  //const int size = n_rows ;
  sicl_gscsrmv_seq(n_rows, vals, c_idx, r_idx, X, R);  // R = A*X;

  // R = B - R
  KERNEL_replace_with_residual_seq(n_rows, B, R) ;
} ;


//
//  S = R + alpha * V
//

static void
smul_vadd_seq(
                  int     n_rows ,
            const float * R      ,
            const float * V      ,
                  float   alpha  ,
                  float * S
         )
{
  //const int size = n_rows;
  KERNEL_smul_vadd_seq( n_rows, R, V, alpha, S) ;
} ;

//
//  X = X + alpha * P
//

static void
add_to_x_seq(
                int     n_rows ,
                float   alpha  ,
          const float * P      ,
                float * X
        )
{
  //const int size = n_rows ;
  KERNEL_add_to_x_seq (n_rows, alpha, P, X);
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
sicl_gscsrmv_seq(
                int     n_rows ,  // numer of matrix rows
          const float * vals   ,  // matrix nonzeros
          const int   * c_idx  ,  // column indices fo nonzeros
          const int   * r_idx  ,  // positions of rows begin, start from 0
          const float * X      ,  // vector, by which matrix is multiplied
                float * R         // result
            )
{
  //const int size = n_rows;
  KERNEL_crs_multiply_seq (
                                              n_rows ,
                                              vals   ,
                                              c_idx  ,
                                              r_idx  ,
                                              X      ,
                                              R
                                            ) ;

  return 0 ;
} ;


int
sicl_gscsrcg_seq(
                  int            n_rows  ,
            const float        * vals    ,
            const int          * c_idx   ,
            const int          * r_idx   ,
                  float        * X       ,
            const float        * B       ,
                  PRECOND_TYPE   precond ,
                  int          * n_iter  ,
                  float        * epsilon
             ){
  int result = -1 ;
  const float almost_zero = numeric_limits<float>::min();

  // Notation, algorithm: see Barlett et al, p. 13
  float rho_1;      // \rho_{i-1}
  float rho_2 = 0;  // \rho_{i-2}
  float alpha = 0;  // \alpha{i}
  float beta;       // \beta_{i-1}

  float norm_b = norm_seq(B, n_rows) ;

  if (norm_b < almost_zero) // if B == 0
  {
    norm_b = 1.0 ;
  }

  float *R, *P, *Z, *Q ;

  // residuals
  R = new float[n_rows]; //cudaCall(cudaMalloc ((void**)(&R), n_rows*sizeof(float)), "cudaMalloc failed for R") ;
  // search directions
  P = new float[n_rows]; //cudaCall(cudaMalloc ((void**)(&P), n_rows*sizeof(float)), "cudaMalloc failed for P") ;
  // solution of A * Z = R
  Z = new float[n_rows]; //cudaCall(cudaMalloc ((void**)(&Z), n_rows*sizeof(float)), "cudaMalloc failed for Z") ;
  Q = new float[n_rows]; //cudaCall(cudaMalloc ((void**)(&Q), n_rows*sizeof(float)), "cudaMalloc failed for Q") ;

  float* pZ = R;
  //pZ = new float[n_rows]; //cudaCall(cudaMalloc ((void**)(&pZ), n_rows*sizeof(float)), "cudaMalloc failed for pZ") ;
  //for (int i=0; i<n_rows; i++)	pZ[i] = R[i]; //cudaCall(cudaMemcpy(pZ, R, n_rows*sizeof(float), cudaMemcpyDeviceToDevice), "cudaMemcpy pZ = R FAILED" ) ;

  // R = B - AX;  (Line 1 at Barlett's algorithm)
  calc_residual_seq(n_rows, vals, c_idx, r_idx, B, X, R) ;

  float residuum = norm_seq(R, n_rows) / norm_b ;

  if (residuum < *epsilon)  // if the trial solution satisfies the equation...
  {
    *n_iter = 0 ;
    result  = 0 ;
  }
  else
  {
    for (int niter = 1 ; niter <= *n_iter ; niter++)
    {
      rho_1 = dot_product_seq(n_rows, R, pZ) ;

      if (1 == niter) {
        // p^1 = z^0;   Barlett: line 6
        for (int i=0; i<n_rows; i++)	P[i] = pZ[i]; //cudaCall(cudaMemcpy(P, pZ, n_rows*sizeof(float), cudaMemcpyDeviceToDevice),          "cudaMemcpy P = pz FAILED" ) ;
      } else {
        beta = (rho_1/rho_2) ;
        smul_vadd_seq(n_rows, pZ, P, beta, P) ;             // P = Z + beta * P
      } ;

      sicl_gscsrmv_seq(n_rows, vals, c_idx, r_idx, P, Q) ; // Q = A*P
      alpha = rho_1 / dot_product_seq(n_rows, P, Q)  ;

      add_to_x_seq(n_rows,  alpha, P, X) ;         // X^i = X^{i-1} + alpha_i * p^i
      add_to_x_seq(n_rows, -alpha, Q, R) ;
      rho_2 = rho_1 ;

      residuum = norm_seq(R, n_rows) / norm_b ;

      if (residuum < *epsilon) // iteration succeeded
      {
        *n_iter = niter ;
        result  = 0 ;
        break ;
      } ;
    } ;
  } ;

  *epsilon = residuum ;

  delete [] Q; //cudaCall(cudaFree (Q), "cudaFree failed for Q") ;
  delete [] R; //cudaCall(cudaFree (R), "cudaFree failed for R") ;
  delete [] P; //cudaCall(cudaFree (P), "cudaFree failed for P") ;
  delete [] Z; //cudaCall(cudaFree (Z), "cudaFree failed for Z") ;
  //delete [] pZ; //cudaCall(cudaFree (pZ), "cudaFree failed for pZ") ;

  return result ;
};

