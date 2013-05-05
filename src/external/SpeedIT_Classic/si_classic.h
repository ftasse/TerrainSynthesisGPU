/*
 * si_classic.h - part of the SpeedIT Classic toolkit
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

#ifndef __SI_CLASSIC_H__
#define __SI_CLASSIC_H__

#ifdef _WIN32

  #ifdef WINDLL
    //#error You should not compile this file as a part of dll library!
    #define DLLAPI __declspec(dllexport)
  #else
    #define DLLAPI //__declspec(dllimport)
  #endif

#else
 #define DLLAPI
#endif

#ifdef __cplusplus
extern "C" {
#endif


////////////////////////////////////////////////////////////////////////////////
//
//  Types of preconditioner used in solver functions (for future expansions)
//
////////////////////////////////////////////////////////////////////////////////
typedef enum {
  P_NONE = 0,
  P_DIAG,
  P_DILU,
  P_DIC
} PRECOND_TYPE ;

////////////////////////////////////////////////////////////////////////////////
//
//  Function solves with preconditioned Conjugate Gradient (CG) method system of
//  linear algebraic equations
//
//      A * X = B
//
//  where A is a sparse matrix and X and B are dense vectors.
//
//
//  Arrays vals, c_idx, r_idx, X and B must be allocated in GPU memory.
//
//  Pointers n_iter and epsilon must point to numbers in CPU memory. These
//  numbers are modified in function.
//
//  ARGUMENTS:
//
//  n_rows  : number of rows and columns in matrix A and positions in 
//            vectors X, R.
//
//  vals    : nonzero elements of matix A. Size of array vals is equal to
//            number of nonzero elements in matrix A.
//
//  c_idx   : column indices for each nonzero element of matrix A. Size of
//            array c_idx is equal to number of nonzero elements of matrix A.
//            Columns are indexed from 0 to n_rows-1.
//
//  r_idx   : indices of row begins for matrix A. Array r_idx contains n_rows+1
//            positions. First position is 0, last position is number of
//            nonzero elements in matrix A.
//
//  X        : dense vectors containing n_rows elements. On input vector must
//             contain trial solution, on output vector is filled with computed
//             solution.
//
//  B        : dense vectors containing n_rows elements. On input vector must 
//             contain right hand side of system of equations.
//
//  precond : type of preconditioner used in solver. 
//            WARNING: preconditioner is UNSUPPORTED in current library 
//            version. Argument declared to keep AP compatibility with
//            future library expansions.
//
//  n_iter   : pointer to maximum number of solver iterations. On input must
//             contain number of maximum allowed iterations, on output is 
//             filled with actually done number of iterations.
//
//  epsilon  : pointer to number containing required solution error. On input
//             must point to numer containing required solution error, on
//             output target is filled with reached solution error.
//
////////////////////////////////////////////////////////////////////////////////
DLLAPI
int 
sicl_gscsrcg(     
                  int            n_rows  ,
            const float        * vals    ,
            const int          * c_idx   ,
            const int          * r_idx   ,
                  float        * X       ,
            const float        * B       ,
                  PRECOND_TYPE   precond ,
                  int          * n_iter  ,
                  float        * epsilon
             ) ;
DLLAPI
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
             ) ;
DLLAPI
int bicgstab_kernel(     int            n_rows  , 
				    const float        * vals    , 
				    const int          * c_idx   , 
				    const int          * r_idx   , 
					  float        * X       , 
				    const float        * B       , 
					  PRECOND_TYPE   precond , // Only for compatibility 
					  int          * n_iter  , 
					  float        * epsilon
             			);
////////////////////////////////////////////////////////////////////////////////
//
//  Function computes sparse matrix by dense vector multiplication
//
//            A * X = R
//
//  where A - square matrix, X and R - vectors. Computation is done with GPU.
//
//  Arrays vals, c_idx, r_idx, X and R must be allocated in GPU memory.
//
//  ARGUMENTS:
//
//  n_rows : number of rows and columns in matrix A and positions in 
//           vectors X, R.
//
//  vals   : nonzero elements of matix A. Size of array vals is equal to
//           number of nonzero elements in matrix A.
//
//  c_idx  : column indices for each nonzero element of matrix A. Size of array
//           c_idx is equal to number of nonzero elements of matrix A. Columns
//           are indexed from 0 to n_rows-1.
//
//  r_idx  : indices of row begins for matrix A. Array r_idx contains n_rows+1
//           positions. First position is 0, last position is number of nonzero
//           elements in matrix A.
//
//  X, R   : dense vectors containing n_rows elements.
//
////////////////////////////////////////////////////////////////////////////////
DLLAPI
int 
sicl_gscsrmv(       
                int     n_rows ,  // numer of matrix rows
          const float * vals   ,  // matrix nonzeros
          const int   * c_idx  ,  // column indices fo nonzeros            
          const int   * r_idx  ,  // positions of rows begin, start from 0 
          const float * X      ,  // vector, by which matrix is multiplied
                float * R         // place for result 
            ) ;
DLLAPI
int 
sicl_gscsrmv_seq(       
                int     n_rows ,  // numer of matrix rows
          const float * vals   ,  // matrix nonzeros
          const int   * c_idx  ,  // column indices fo nonzeros            
          const int   * r_idx  ,  // positions of rows begin, start from 0 
          const float * X      ,  // vector, by which matrix is multiplied
                float * R         // place for result 
            ) ;            

#ifdef __cplusplus
} ;
#endif


#endif
