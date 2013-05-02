#include "dtts_2dmorph.h"

TPS_Morpher::TPS_Morpher(
    std::auto_ptr<std::vector<Coord_Diff> > p_samples,
    float regularization,
    float subsampling_factor)
    : samples( p_samples )
{
    assert( samples->size() >= 3 );
    unsigned p = samples->size();

    unsigned m = (unsigned)(p * subsampling_factor);
    if ( m < 3 ) m = 3;
    if ( m > p ) m = p;

    if ( m < p )
    {
        // Randomize the input if subsampling is used
        for ( unsigned i=0; i<samples->size()-1; ++i )
        {
            int j = i + ((unsigned)rand()) % (samples->size()-i);
            Coord_Diff tmp = (*samples)[j];
            (*samples)[i] = (*samples)[j];
            (*samples)[j] = tmp;
        }
    }

    // Allocate the matrix and vector
    mtx_l.resize(p+3, m+3);
    mtx_v.resize(p+3, 2);
    mtx_orig_k.resize(p, m);

    // Fill K (p x m, upper left of L)
    for ( unsigned i=0; i<p; ++i )
    {
        const Coord_Diff& pt_i = (*samples)[i];
        for ( unsigned j=0; j<m; ++j )
        {
            const Coord_Diff& pt_j = (*samples)[j];
            float elen2 = SQR(pt_i.x-pt_j.x) + SQR(pt_i.y-pt_j.y);
            mtx_l(i,j) = mtx_orig_k(i,j) = base_func(elen2);
        }
    }

    // Empiric value for avg. distance between points
    //
    // This variable is normally calculated to make regularization
    // scale independent, but since our shapes in this application are always
    // normalized to maxspect [-.5,.5]x[-.5,.5], this approximation is pretty
    // safe and saves us p*p square roots
    const float a = 0.5;

    // Fill the rest of L
    for ( unsigned i=0; i<p; ++i )
    {
        const Coord_Diff pt_i = (*samples)[i];

        // P (p x 3, upper right)
        mtx_l(i, m+0) = 1.0;
        mtx_l(i, m+1) = pt_i.x;
        mtx_l(i, m+2) = pt_i.y;

        if ( i<m )
        {
            // diagonal: reqularization parameters (lambda * a^2)
            mtx_l(i,i) = mtx_orig_k(i,i) =
                             regularization * (a*a);

            // P transposed (3 x p, bottom left)
            mtx_l(p+0, i) = 1.0;
            mtx_l(p+1, i) = pt_i.x;
            mtx_l(p+2, i) = pt_i.y;
        }
    }

    // O (3 x 3, lower right)
    for ( unsigned i=p; i<p+3; ++i )
        for ( unsigned j=m; j<m+3; ++j )
            mtx_l(i,j) = 0.0;

    // Fill the right hand matrix V
    for ( unsigned i=0; i<p; ++i )
    {
        const Coord_Diff& pt_i = (*samples)[i];
        mtx_v(i,0) = pt_i.dx;
        mtx_v(i,1) = pt_i.dy;
    }

    mtx_v(p+0, 0) = mtx_v(p+1, 0) = mtx_v(p+2, 0) = 0.0;
    mtx_v(p+0, 1) = mtx_v(p+1, 1) = mtx_v(p+2, 1) = 0.0;

    // Solve the linear system "inplace"
    int sret = LU_Solve(mtx_l, mtx_v);
    assert( sret != 2 );
    if (sret == 1)
    {
        //puts( "Singular matrix! Aborting." );
        //exit(1);
        singularmatrix=true;
    }
    else singularmatrix=false;
}

void TPS_Morpher::morph( std::vector<Point>& pts )
{
    //assert( samples.get() && "Morpher no longer owns 'samples'");

    const unsigned m = mtx_orig_k.size2();
    for ( std::vector<Point>::iterator ite=pts.begin(), end=pts.end();
            ite != end;
            ++ite )
    {
        float x=ite->x, y=ite->y;
        float dx = mtx_v(m+0, 0) + mtx_v(m+1, 0)*x + mtx_v(m+2, 0)*y;
        float dy = mtx_v(m+0, 1) + mtx_v(m+1, 1)*x + mtx_v(m+2, 1)*y;

        std::vector<Coord_Diff>::const_iterator diff_ite = samples->begin();
        Matrix_Col cv0(mtx_v,0), cv1(mtx_v,1);
        Matrix_Col::const_iterator cv0_ite(cv0.begin()), cv1_ite(cv1.begin());
        for ( unsigned i=0; i<m; ++i, ++diff_ite, ++cv0_ite, ++cv1_ite )
        {
            float d = base_func( SQR(diff_ite->x - x) + SQR(diff_ite->y - y) );
            dx += (*cv0_ite) * d;
            dy += (*cv1_ite) * d;
        }

        ite->x += dx;
        ite->y += dy;
    }
}

void TPS_Morpher::morph_pt( Point& pts )
{
    //assert( samples.get() && "Morpher no longer owns 'samples'");

    const unsigned m = mtx_orig_k.size2();
    float x=pts.x, y=pts.y;
    float dx = mtx_v(m+0, 0) + mtx_v(m+1, 0)*x + mtx_v(m+2, 0)*y;
    float dy = mtx_v(m+0, 1) + mtx_v(m+1, 1)*x + mtx_v(m+2, 1)*y;

    std::vector<Coord_Diff>::const_iterator diff_ite = samples->begin();
    Matrix_Col cv0(mtx_v,0), cv1(mtx_v,1);
    Matrix_Col::const_iterator cv0_ite(cv0.begin()), cv1_ite(cv1.begin());
    for ( unsigned i=0; i<m; ++i, ++diff_ite, ++cv0_ite, ++cv1_ite )
    {
        float d = base_func( SQR(diff_ite->x - x) + SQR(diff_ite->y - y) );
        dx += (*cv0_ite) * d;
        dy += (*cv1_ite) * d;
    }

    pts.x += dx;
    pts.y += dy;
}

float TPS_Morpher::calc_bending_energy()
{
    //assert( samples.get() && "Morpher no longer owns 'samples'");

    // bending energy = trace( W^T * A * W ),
    // where A = upper left m x m block of mtx_orig_k
    const unsigned m = mtx_orig_k.size2();
    ublas::matrix_range<Matrix> mtx_w(mtx_v,
                                      ublas::range(0, m),
                                      ublas::range(0, 2));
    ublas::matrix_range<Matrix> mtx_a(mtx_orig_k,
                                      ublas::range(0, m),
                                      ublas::range(0, m));
    Matrix bm = prod( Matrix(prod(trans(mtx_w), mtx_a)), mtx_w);
    assert( bm.size1() == bm.size2() && bm.size1() == 2 );
    return bm(0,0) + bm(1,1);
}

std::auto_ptr<std::vector<Coord_Diff> > TPS_Morpher::grab_samples()
{
    return samples;
}

inline float TPS_Morpher::base_func(float r2)
{
    // same as r*r * logf(r), but for r^2:
    return ( r2==0 )
           ? 0.0 // function limit at 0
           : r2 * logf(r2) * 0.217147241; // = 1/(2*logf(10))
}
