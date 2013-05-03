#ifndef DTTS_2DMORPH_H
#define DTTS_2DMORPH_H

/**
 *  Thin Plate Spline 2D point morpher example.
 *
 *    Takes in sample point coordinates and their displacements,
 *    fits a TPS into them and allowes morphing new points
 *    accordint to the approximated deformation.
 *
 *    Supports TPS approximation 3 as suggested in paper
 *    Gianluca Donato and Serge Belongie, 2002: "Approximation
 *    Methods for Thin Plate Spline Mappings and Principal Warps"
 *
 *  This should be considered more as an example than a ready made module!
 *  The code has been extracted from a working "deformed shape matching"
 *  application and has been optimized for that particular case.
 *  I don't even know if this compiles or not.
 *
 *  Copyright (C) 2003-2005 by Jarno Elonen
 *
 *  This is Free Software / Open Source with a very permissive
 *  license:
 *
 *  Permission to use, copy, modify, distribute and sell this software
 *  and its documentation for any purpose is hereby granted without fee,
 *  provided that the above copyright notice appear in all copies and
 *  that both that copyright notice and this permission notice appear
 *  in supporting documentation.  The authors make no representations
 *  about the suitability of this software for any purpose.
 *  It is provided "as is" without express or implied warranty.
 */

#include "external/ludecomposition.h"

#include <vector>
#include <utility>

#include <cassert>
#include <cmath>

using namespace boost::numeric;
typedef ublas::matrix<float> Matrix;
typedef ublas::matrix_row<Matrix> Matrix_Row;
typedef ublas::matrix_column<Matrix> Matrix_Col;


inline float SQR( float x )
{
    return x*x;
}

/// A 2D point
struct Point
{
    float x, y;
};

/// A displacement of one 2D point: x,y is the original location
/// and dx,dy is the offset.
struct Coord_Diff
{
    Coord_Diff()
    {}

    Coord_Diff( float x, float y, float dx, float dy )
        : x(x), y(y), dx(dx), dy(dy)
    {}

    float x, y, dx, dy;
};

/// 2D point morph interpolator using TPS (Thin Plate Spline).
class TPS_Morpher
{
public:

    /// Calculate the morph weights from sample points.
    /// Builds a matrix of (approximately) size N(p_samples) x N(p_samples*subsampling_factor),
    /// and inverts it using LU decomposition (O(N^3), so be careful not to input too large sample sets.
    ///
    /// For performance reasons, this function assumes ownership of the sample vector so don't
    /// change it after passing it here. Once you have done using the morpher, you can claim
    /// it back by calling grab_samples().
    ///
    /// The function assumes that the average distance between the points to be about 0.5 to
    /// save some computation. If this is not the case in your application, you need to calculate
    /// (or approximate) this value by yourself. See the comments concerning "float a" in the code.
    ///
    /// @param p_samples Displacement samples to be interpolated
    /// @param regularization Amount of "relaxation", 0.0 = exact interpolation
    /// @param subsampling_factor 1.0 = use all points, 0.5 = use 1/2 of the points etc.
    TPS_Morpher(
        std::auto_ptr<std::vector<Coord_Diff> > p_samples,
        float regularization,
        float subsampling_factor = 1.0 );

    /// Morph given points according to the TPS
    /// @param pts The points to morph
    void morph_pt( Point& pts );
    void morph( std::vector<Point>& pts );

    /// Calculate bending energy, or if subsampling_factor
    /// for constructor was < 1.0, approximate it.
    float calc_bending_energy();

    /// Takes away the ownership of 'samples' from the morpher.
    /// After this, calling other functions becomes illegal.
    std::auto_ptr<std::vector<Coord_Diff> > grab_samples();

    bool singularmatrix;
private:

    inline float base_func(float r2);

    std::auto_ptr<std::vector<Coord_Diff> > samples;
    Matrix mtx_l, mtx_v, mtx_orig_k;
};

#endif
