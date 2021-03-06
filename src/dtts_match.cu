#include "cutil_math_bugfixes.h"
#include "dtts_match.incl"

#include <thrust/sort.h>
#include <thrust/device_ptr.h>
#include <thrust/device_malloc.h>
#include <thrust/device_free.h>
#include <thrust/device_vector.h>
#include <thrust/system_error.h>
#include <thrust/count.h>

#define BLOCK_SIZE 256

texture<float, 2, cudaReadModeElementType> src_tex;
texture<float, 2, cudaReadModeElementType> utar_tex;
extern __device__ int bigw = 0;

#if defined (__CUDACC__)
__host__ __device__ 
#endif
inline float2 operator-(float a, float2 b)
{
    return make_float2(a - b.x, a - b.y);
}

inline __device__ void bspline_weights(float2 fraction, float2& w0, float2& w1, float2& w2, float2& w3)
{
    const float2 one_frac = 1.0f-fraction;
    const float2 squared = fraction * fraction;
    const float2 one_sqd = one_frac * one_frac;

    w0 = 1.0f/6.0f * one_sqd * one_frac;
    w1 = 2.0f/3.0f-0.5f * squared * (2.0f-fraction);
    w2 = 2.0f/3.0f-0.5f * one_sqd * (2.0f-one_frac);
    w3 = 1.0f/6.0f * squared * fraction;
}

__device__ float cubicInterpol(int mwidth, int mheight, const float fx, const float fy)
{
    // transform the coordinate from [0,extent] to [-0.5, extent-0.5]
    /*const float2 coord_grid = make_float2(fx - 0.5, fy - 0.5);
    float2 index = floor(coord_grid);
    const float2 fraction = coord_grid - index;
    index.x += 0.5;  //move from [-0.5, extent-0.5] to [0, extent]
    index.y += 0.5;  //move from [-0.5, extent-0.5] to [0, extent]

    float result = 0.0;
    for (float y=-1; y < 2.5; y++)
    {
    	float bsplineY = bspline(y-fraction.y);
    	float v = index.y + y;
    	for (float x=-1; x < 2.5; x++)
    	{
    		float bsplineXY = bspline(x-fraction.x) * bsplineY;
    		float u = index.x + x;
    		result += bsplineXY * tex2D(src_tex, u, v);
    	}
    }
    return result;*/

    /*const float2 coord_grid = make_float2(fx - 0.5f, fy - 0.5f);
    const float2 index = floor(coord_grid);
    const float2 fraction = coord_grid - index;
    float2 w0, w1, w2, w3;
    bspline_weights(fraction, w0, w1, w2, w3);

    const float2 g0 = w0 + w1;
    const float2 g1 = w2 + w3;
    const float2 h0 = (w1 / g0) - make_float2(0.5f) + index;  //h0 = w1/g0 - 1, move from [-0.5, extent-0.5] to [0, extent]
    const float2 h1 = (w3 / g1) + make_float2(1.5f) + index;  //h1 = w3/g1 + 1, move from [-0.5, extent-0.5] to [0, extent]

    // fetch the four linear interpolations
    float tex00 = tex2D(src_tex, h0.x, h0.y);
    float tex10 = tex2D(src_tex, h1.x, h0.y);
    float tex01 = tex2D(src_tex, h0.x, h1.y);
    float tex11 = tex2D(src_tex, h1.x, h1.y);
    // weigh along the y-direction
    tex00 = g0.y * tex00 + g1.y * tex01;
    tex10 = g0.y * tex10 + g1.y * tex11;

    // weigh along the x-direction
    return (g0.x * tex00 + g1.x * tex10);*/

    const float
    nfx = fx,//<0?0:(fx>mwidth-1?mwidth-1:fx),
          nfy = fy;//<0?0:(fy>mheight-1?mheight-1:fy);
    const int
    x = (int) nfx,
        y = (int) nfy;
    const float
    dx = nfx-x,
         dy = nfy-y;

    const int
    px = x-1, //<0?0:x-1,
         nx = x+1, //dx>0?x+1:x,
              ax = x+2, //>=mwidth?mwidth-1:x+2,
                   py = y-1, //<0?0:y-1,
                        ny = y+1,//dy>0?y+1:y,
                             ay = y+2; //>=mheight?mheight-1:y+2;

    const float
    Ipp = tex2D(src_tex,px,py), Icp = tex2D(src_tex,x,py), Inp = tex2D(src_tex,nx,py), Iap = tex2D(src_tex,ax,py),
                                      Ip = Icp + 0.5f*(dx*(-Ipp+Inp) + dx*dx*(2*Ipp-5*Icp+4*Inp-Iap) + dx*dx*dx*(-Ipp+3*Icp-3*Inp+Iap)),
                                           Ipc = tex2D(src_tex,px,y),  Icc = tex2D(src_tex,x,y), Inc = tex2D(src_tex,nx,y),  Iac = tex2D(src_tex,ax,y),
                                                   Ic = Icc + 0.5f*(dx*(-Ipc+Inc) + dx*dx*(2*Ipc-5*Icc+4*Inc-Iac) + dx*dx*dx*(-Ipc+3*Icc-3*Inc+Iac)),
                                                           Ipn = tex2D(src_tex,px,ny), Icn = tex2D(src_tex,x,ny), Inn = tex2D(src_tex,nx,ny), Ian = tex2D(src_tex,ax,ny),
                                                                   In = Icn + 0.5f*(dx*(-Ipn+Inn) + dx*dx*(2*Ipn-5*Icn+4*Inn-Ian) + dx*dx*dx*(-Ipn+3*Icn-3*Inn+Ian)),
                                                                           Ipa = tex2D(src_tex,px,ay), Ica = tex2D(src_tex,x,ay), Ina = tex2D(src_tex,nx,ay), Iaa = tex2D(src_tex,ax,ay),
                                                                                   Ia = Ica + 0.5f*(dx*(-Ipa+Ina) + dx*dx*(2*Ipa-5*Ica+4*Ina-Iaa) + dx*dx*dx*(-Ipa+3*Ica-3*Ina+Iaa));

    return Ic + 0.5f*(dy*(-Ip+In) + dy*dy*(2*Ip-5*Ic+4*In-Ia) + dy*dy*dy*(-Ip+3*Ic-3*In+Ia));
}

__device__ float getH_checkBounds(register int src_w, register int src_h, register int bsize, register int rx, register int ry, register int rt, register int i, register int j )
{
    int ri = 0;
    int rj = 0;
    float candv=0;
    int mid = bsize/2;

    if (rt<360/DROT)
    {
        float ni=0;
        float nj=0;
        int rot = rt*DROT;
        ni = (rx+mid) + ((i - mid)*cosf(rot*180.0f/PI)) + ((j - mid) *sinf(rot*180.0f/PI));
        nj = (ry+mid) - ((i - mid)*sinf(rot*180.0f/PI)) + ((j - mid) *cosf(rot*180.0f/PI));
        if (ni>=0 && nj>=0 && ni<src_w && nj<src_h)
            candv = cubicInterpol(src_w,src_h,ni,nj);
        return candv;
    }
    else if (rt==360/DROT)
    {
        ri = (rx+bsize-i-1);
        rj = (ry+j);
        if (ri>=0 && rj>=0 && ri<src_w && rj<src_h)
            return  tex2D(src_tex,ri , rj );
    }
    else
    {
        ri = (rx+i);
        rj = (ry+bsize-j-1);
        if (ri>=0 && rj>=0 && ri<src_w && rj<src_h)
            return  tex2D(src_tex,ri , rj );
    }

    return candv;
}

__device__ float getH(int src_w, int src_h, int bsize, int rx, int ry, int rt, int i, int j )
{
    int ri = 0;
    int rj = 0;
    float candv=0;
    int mid = bsize/2;

    if (rt<360/DROT)
    {
        float ni=0;
        float nj=0;
        int rot = rt*DROT;
        ni = (rx+mid) + ((i - mid)*cosf(rot*180.0f/PI)) + ((j - mid) *sinf(rot*180.0f/PI));
        nj = (ry+mid) - ((i - mid)*sinf(rot*180.0f/PI)) + ((j - mid) *cosf(rot*180.0f/PI));
        candv = cubicInterpol(src_w,src_h,ni,nj);
        return candv;
    }
    else if (rt==360/DROT)
    {
        ri = (rx+bsize-i-1);
        rj = (ry+j);
        return  tex2D(src_tex,ri , rj );
    }
    else
    {
        ri = (rx+i);
        rj = (ry+bsize-j-1);
        return  tex2D(src_tex,ri , rj );
    }

    //return candv;
}

__device__ float get_diff(int src_w, int src_h, int bsize,  node_t pv, int rt, point_t p, point_t q, int dx, int dy)
{
    int px = pv.x;
    int py = pv.y;
    //int alln = bigwidth/(bsize*bsize);
    float sum=0;
    int count=0;
    int x0=p.x, y0=p.y;
    int x1=q.x, y1=q.y;

    bool steep = ( abs(y1 - y0) > abs(x1 - x0) );
    if (steep)
    {
        swap(x0, y0);
        swap(x1, y1);
    }
    if (x0 > x1)
    {
        swap(x0, x1);
        swap(y0, y1);
    }
    int deltax = x1 - x0;
    int deltay = abs(y1 - y0);
    int error = deltax / 2 ;
    int ystep;
    int y = y0;
    if (y0 < y1)
        ystep = 1 ;
    else ystep = -1;

    for (int x=x0; x<=x1; x++)
    {
        //if (x>=0 && x<bsize && y>=0 && y<bsize )
        {
            if (steep)
            {
                //int ide = y + x*bsize;
                float sval = getH(src_w,src_h, bsize, px,py,rt,y,x);
                float val = tex2D(utar_tex, y+dx, x+dy)-sval;
                sum+= val*val;
            }
            else
            {
                //int ide = x + y*bsize;
                float sval = getH(src_w,src_h, bsize, px,py,rt,x,y);
                float val = tex2D(utar_tex, x+dx, y+dy)-sval;

                sum+= val*val;
            }
            count++;
            error = error - deltay;
            if (error < 0)
            {
                y = y + ystep;
                error = error + deltax;
            }
        }
    }

    return sqrtf(sum)/count;
}

__device__ float get_diff_profile(int sw, int sh, int bsize, node_t pv, int rt, node_t* cand_br1, int csize, int dx, int dy)
{
    /*float total=0.;
    int mw=bsize, mh=bsize;
    node_t mid(mw/2,mh/2);

    for (int k=0; k<csize; k++)
    {
        total += get_diff(sw,sh,bsize,pv,rt,cand_br1[k],mid,dx,dy);
    }
    if (csize>2)  return total;

    if (csize==2)
    {
        node_t qnode1,qnode2;

        qnode1 = cand_br1[0];
        qnode2=cand_br1[1];
        node_t vec(-(qnode1.y-qnode2.y)/2,(qnode1.x-qnode2.x)/2);

        qnode1 = node_t(mid.x+vec.x,mid.y+vec.y);
        qnode2 = node_t(mid.x-vec.x,mid.y-vec.y);
        total += get_diff(sw,sh,bsize,pv,rt,qnode1, mid,dx,dy);
        total += get_diff(sw,sh,bsize,pv,rt,mid, qnode2,dx,dy);
    }

    return total;*/

    float total=0.;
    int mw=bsize, mh=bsize;
    point_t mid(mw/2,mh/2);


    //if (csize>2)  return total;

    //if (csize==2)
    {

        for (int k=0; k<csize; k++){

            point_t qnode1;
            qnode1 = point_t(cand_br1[k].x, cand_br1[k].y) ;
            point_t vec1(qnode1.x-mid.x, qnode1.y-mid.y);

            point_t pv1(-vec1.y, vec1.x);

            point_t p1, np1;

            p1 = point_t(mid.x+pv1.x, mid.y+pv1.y);
            np1 = point_t(mid.x-pv1.x, mid.y-pv1.y);
            total += get_diff(sw,sh,bsize,pv,rt,p1,np1, dx, dy);

            p1 = point_t(qnode1.x+pv1.x, qnode1.y+pv1.y);
            np1 = point_t(qnode1.x-pv1.x, qnode1.y-pv1.y);
            total += get_diff(sw,sh,bsize,pv,rt,p1,np1, dx, dy);

            qnode1.x =((mid.x+cand_br1[k].x)*1)/2; qnode1.y =((mid.y+cand_br1[k].y)*1)/2;

            p1 = point_t(qnode1.x+pv1.x, qnode1.y+pv1.y);
            np1 = point_t(qnode1.x-pv1.x, qnode1.y-pv1.y);
            total += get_diff(sw,sh,bsize,pv,rt,p1,np1, dx, dy);
        }

        //total/=3;
        //total = 0;

         for (int k=0; k<csize; k++)
        {
            point_t qnode1;
            qnode1 = point_t(cand_br1[k].x, cand_br1[k].y) ;
            total += get_diff(sw,sh,bsize,pv,rt,qnode1,mid, dx, dy);
        }
    }
        return total/csize;
}


__device__ float get_diff_profile(int sw, int sh, int bsize, node_t pv, int rt, node_t* cand_br1, int csize, int dx, int dy)
{
    /*float total=0.;
    int mw=bsize, mh=bsize;
    node_t mid(mw/2,mh/2);

    for (int k=0; k<csize; k++)
    {
        total += get_diff(sw,sh,bsize,pv,rt,cand_br1[k],mid,dx,dy);
    }
    if (csize>2)  return total;

    if (csize==2)
    {
        node_t qnode1,qnode2;

        qnode1 = cand_br1[0];
        qnode2=cand_br1[1];
        node_t vec(-(qnode1.y-qnode2.y)/2,(qnode1.x-qnode2.x)/2);

        qnode1 = node_t(mid.x+vec.x,mid.y+vec.y);
        qnode2 = node_t(mid.x-vec.x,mid.y-vec.y);
        total += get_diff(sw,sh,bsize,pv,rt,qnode1, mid,dx,dy);
        total += get_diff(sw,sh,bsize,pv,rt,mid, qnode2,dx,dy);
    }

    return total;*/

    float total=0.;
    int mw=bsize, mh=bsize;
    point_t mid(mw/2,mh/2);


    //if (csize>2)  return total;

    //if (csize==2)
    {

        for (int k=0; k<csize; k++){

            point_t qnode1;
            qnode1 = point_t(cand_br1[k].x, cand_br1[k].y) ;
            point_t vec1(qnode1.x-mid.x, qnode1.y-mid.y);

            point_t pv1(-vec1.y, vec1.x);

            point_t p1, np1;

            p1 = point_t(mid.x+pv1.x, mid.y+pv1.y);
            np1 = point_t(mid.x-pv1.x, mid.y-pv1.y);
            total += get_diff(sw,sh,bsize,pv,rt,p1,np1, dx, dy);

            p1 = point_t(qnode1.x+pv1.x, qnode1.y+pv1.y);
            np1 = point_t(qnode1.x-pv1.x, qnode1.y-pv1.y);
            total += get_diff(sw,sh,bsize,pv,rt,p1,np1, dx, dy);

            qnode1.x =((mid.x+cand_br1[k].x)*1)/2; qnode1.y =((mid.y+cand_br1[k].y)*1)/2;

            p1 = point_t(qnode1.x+pv1.x, qnode1.y+pv1.y);
            np1 = point_t(qnode1.x-pv1.x, qnode1.y-pv1.y);
            total += get_diff(sw,sh,bsize,pv,rt,p1,np1, dx, dy);
        }

        //total/=3;
        //total = 0;

         for (int k=0; k<csize; k++)
        {
            point_t qnode1;
            qnode1 = point_t(cand_br1[k].x, cand_br1[k].y) ;
            total += get_diff(sw,sh,bsize,pv,rt,qnode1,mid, dx, dy);
        }
    }
        return total/csize;
}

__device__ float ssdf(float* dest, int src_w, int src_h, int bsize, int px, int py, int rt)
{
    float sum = 0.;
    int count = 0;
    for (int x=0; x<bsize*bsize; x++)
    {
        float dval = dest[x];
        if (dval!=0.)
        {
            float sval = getH(src_w,src_h, bsize, px,py,rt,(x%bsize),(x/bsize));
            float val = dval-sval;
            sum+= val*val;
            count++;
        }
    }
    count++;
    return sqrtf(sum)/count;
}

__global__ void fillCandidates_kernel(float* cands, cost_t* choices, int num, int bsize, int src_w, int src_h )
{
    unsigned int tid = blockDim.x * blockIdx.x + threadIdx.x;
    if (tid<num*bsize*bsize)
    {
        int k = (tid)/(bsize*bsize);
        int id = (tid)%(bsize*bsize);
        cands[tid] = getH(src_w,src_h,bsize,choices[k].org.x,choices[k].org.y,choices[k].rot/DROT,id%bsize,id/bsize);
    }
}

Image findCand(Image& dest, vector<cost_t> candidates, int src_w, int src_h, cost_t& prev, int bsize, int dx, int dy)
{
    //cout<<dx<<" "<<dy<<endl;
    int s = candidates.size();
    if (s>KNUM) s=KNUM;

    vector<cost_t> choices;
    //cout<<prev.org.x<<" "<<prev.org.y<<" "<<prev.rot<<" "<<s<<endl;
    int k=0, it=0, csize=candidates.size();
    while (k<csize && it<s)	{
        if (prev.org.x!=candidates[k].org.x || prev.org.y!=candidates[k].org.y )
        {
            //cout<<candidates[k].org.x<<" "<<candidates[k].org.y<<" "<<candidates[k].rot<<endl;
            choices.push_back(candidates[k]);
            k++;
        }
        it++;
    }
    if (choices.size()==0)  choices.push_back(candidates[0]);

    float* cand_vals = new float [choices.size()*bsize*bsize];
    float* cand_vals_dev;
    cudaMalloc((void**) &cand_vals_dev, sizeof(float)*choices.size()*bsize*bsize);
    cost_t* choices_dev;
    cudaMalloc((void**) &choices_dev, sizeof(cost_t)*choices.size());
    cudaMemcpy(choices_dev, &choices[0], sizeof(cost_t)*choices.size(),cudaMemcpyHostToDevice);

    int num_blocks = (choices.size()*bsize*bsize)/BLOCK_SIZE+1;
    fillCandidates_kernel<<<num_blocks,BLOCK_SIZE>>>( cand_vals_dev, choices_dev, choices.size(), bsize, src_w, src_h );
    cudaMemcpy(cand_vals, cand_vals_dev, sizeof(float)*choices.size()*bsize*bsize,cudaMemcpyDeviceToHost);
    cudaFree(cand_vals_dev);
    cudaFree(choices_dev);

    float mini = 10*INF;
    int ch = 0;

    for (unsigned int k=0; k<choices.size(); k++)
    {
        Image cand(bsize,bsize);
        getCand(cand,cand_vals,k);
        float cost = graphCut_cost(&dest,&cand,dx,dy);
        if (cost<mini)
        {
            mini = cost;
            ch = k;
        }
    }

    Image res(bsize,bsize);
    getCand(res,cand_vals,ch);

    prev = choices[ch];
    delete [] cand_vals;

    // cout<<"Choice: "<<choices[ch].org.x<<" "<<choices[ch].org.y<<" "<<choices[ch].rot<<endl;


    return res;
}

__global__ void buildCandidates_noFeature_kernel(cost_t* candidates, node_t* dnodes, int src_w, int src_h, int nsize,  int bsize)
{


    unsigned int k = blockDim.x * blockIdx.x + threadIdx.x;
    int rs = (360/DROT+ DMIR);
    if (k<(nsize*rs))
    {


        int kn = k%(nsize);
        int rt = k/(nsize);
        node_t pnode = dnodes[kn];
        int rx = pnode.x;
        int ry = pnode.y;

        int rot = rt*DROT;

        cost_t c(node_t(rx,ry),0,rot,0,0);

        {
            //int kpos = (kn+rt*nsize);
            c.vpos = kn*3;

            c.skip = false;

            for (int id=0; id<bsize*bsize; id++)
            {

                int i = id%bsize;
                int j = id/bsize;
                //cout<<i<<"/"<i<<": "<<getH(src.getPixels(),src.width(),src.height(),bsize,rx,ry,rt,i,j)<<endl;
                float candv= getH_checkBounds(src_w,src_h,bsize,rx,ry,rt,i,j);
                if (candv<=BG)
                {
                    c.skip = true;
                    //break;
                }
            }

            candidates[k]=c;



        }
    }
}

__device__ float getCost_noFeature(float* dest, float* dem_vars, float* usr_var, cost_t dem, int sw, int sh, int bsize)
{
    node_t p= dem.org;
    int rt = dem.rot/DROT;
    float tmp = 0;
    tmp+=50*ssdf(dest, sw,sh,bsize,p.x,p.y,rt);
    if (use_noisestat) tmp+=0.0001*compare_variances(usr_var,dem_vars,dem.vpos,NLEVEL);
    return tmp;
}

__global__ void ComputeCosts_noFeature_kernel(float* dest, int sw, int sh, cost_t* candidates, int csize, float* dem_vars, float* usr_var, int bsize)
{
    int k = blockDim.x * blockIdx.x + threadIdx.x;
    if (k<csize)
    {
        if ( (!candidates[k].skip))
        {
            candidates[k].cost = getCost_noFeature(dest, dem_vars, usr_var, candidates[k], sw, sh, bsize);
        }
        else
        {
            candidates[k].cost = INF;
        }
    }
}

void match_noFeature_gpu2(Terrain& dest, Image& src, Image& target, node_list dem_nodes, vector<Image>& tar_pyr, vector<Image>&  src_pyr, int bsize, int osize)
{
    clock_t start_t, end_t, s_tmp, e_tmp;
    start_t = clock();

    match_time = 0;
    paste_time = 0;
    get_target = 0;
    find_time = 0;
    sort_time = 0;

    int nsize = dem_nodes.size();
    int rs = 360/DROT;

    //node_t* dem_leafs = new node_t [nsize*(rs+DMIR)*10];
    //node_t* dem_lsizes = new node_t [nsize*(rs+DMIR)];
    node_t* dnodes = new node_t [nsize];


    s_tmp = clock();

    float* dem_vars =  matchingPrepocessing(src_pyr, dem_nodes, dnodes,bsize);

    {
        e_tmp = clock();
        float elapsed = ((float)( e_tmp - s_tmp )) /CLOCKS_PER_SEC;
        cerr<<"Get noise stats and control points elapsed time: "<<elapsed<<" s.\n";
        s_tmp = clock();
    }

    float cand_time = 0;
    //cout<<"Prepare candidates! "<< nsize<<"\n";

    thrust::device_vector<cost_t> candidates_dev(nsize*(rs+DMIR));
    node_t* dnodes_dev;
    cudaMalloc((void**) &dnodes_dev, sizeof(node_t)*nsize);
    cudaMemcpy(dnodes_dev, &dnodes[0], sizeof(node_t)*nsize, cudaMemcpyHostToDevice);

    float* dem_vars_dev;
    cudaMalloc((void**) &dem_vars_dev, sizeof(float)*nsize*NLEVEL);
    cudaMemcpy(dem_vars_dev, dem_vars, sizeof(float)*nsize*NLEVEL, cudaMemcpyHostToDevice);

    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
    cudaArray* src_dev;
    cudaMallocArray(&src_dev, &channelDesc, src.width(), src.height());
    cudaMemcpyToArray(src_dev, 0, 0, src.getPixels(), sizeof(float)*src.width()*src.height(), cudaMemcpyHostToDevice);
    // Set texture parameters
    src_tex.addressMode[0] = cudaAddressModeClamp;
    src_tex.addressMode[1] = cudaAddressModeClamp;
    src_tex.filterMode = cudaFilterModePoint;
    src_tex.normalized = false;
    cudaBindTextureToArray(src_tex, src_dev, channelDesc);

    cudaArray* target_dev;
    cudaMallocArray(&target_dev, &channelDesc, target.width(), target.height());
    cudaMemcpyToArray(target_dev, 0, 0, target.getPixels(), sizeof(float)*target.width()*target.height(), cudaMemcpyHostToDevice);
    // Set texture parameters
    utar_tex.addressMode[0] = cudaAddressModeClamp;
    utar_tex.addressMode[1] = cudaAddressModeClamp;
    utar_tex.filterMode = cudaFilterModePoint;
    utar_tex.normalized = false;
    cudaBindTextureToArray(utar_tex, target_dev, channelDesc);

    int threadsPerBlock = BLOCK_SIZE;
    int blocksPerGrid =  ((nsize*(rs+DMIR)) / threadsPerBlock)+1;
    {
        buildCandidates_noFeature_kernel<<<blocksPerGrid, threadsPerBlock>>>(thrust::raw_pointer_cast(&candidates_dev[0]), &dnodes_dev[0], src.width(), src.height(), nsize, bsize);
        cudaThreadSynchronize();
        e_tmp = clock();
        float elapsed = ((float)( e_tmp - s_tmp )) /CLOCKS_PER_SEC;
        //cerr<<"Build non-feature matching candidates elapsed time: "<<elapsed<<" s.\n";
        cand_time += elapsed;
    }

    cost_t prev (node_t(-1,-1),0,0,0,0);
    int dx = -5*bsize, dy=-5*bsize;

    int cnum = nsize*(rs+DMIR);
    cout<<"Start matching non-feature patches! from "<<cnum<<"\n";

    float* usr_var_dev;
    cudaMalloc((void**) &usr_var_dev,sizeof(float)*NLEVEL);
    float* ucand_dev;
    cudaMalloc((void **) &ucand_dev, sizeof(float)*bsize*bsize);
    cnum = 0;

    vector<cost_t> candidates (nsize*(rs+DMIR));

    list<cost_t> omega;
    float* conf = new float [dest.width()*dest.height()];
    for (int x=0; x<dest.width(); x++)	for (int y=0; y<dest.height(); y++)
            if (dest(x,y)>BG)
                conf[x+y*dest.width()] = 1;
            else
                conf[x+y*dest.width()] = 0;
    while (true)
    {
        s_tmp = clock();

        getNextTarget (dest.getPixels(), conf, omega,dest.width(), dest.height(), bsize, dx+bsize/2, dy+bsize/2);
        if (omega.size()==0) break;


        int xmaxi = -1;
        int ymaxi = -1;
        float maxv = -1e6f;
        float tmp=0.0f;

        for (list<cost_t>::iterator it=omega.begin(); it!=omega.end(); it++)
        {
            tmp = (*it).cost;
            if (tmp>maxv+1e-6f)
            {
                node_t p = (*it).org;
                maxv = tmp;
                xmaxi = p.x;
                ymaxi = p.y;
            }
        }

        e_tmp = clock();
        get_target+=mstimer(s_tmp,e_tmp);

        dx=(xmaxi)-bsize/2;
        dy=(ymaxi)-bsize/2;


        vector<float> usr_var = noise_variances(tar_pyr,dx,dy,bsize);
        Image ucand = dest.get_crop(dx,dy,dx+bsize-1,dy+bsize-1);
        //cout<<"hi0\n";

        s_tmp = clock();
        cudaMemcpy(usr_var_dev, &usr_var[0], sizeof(float)*NLEVEL, cudaMemcpyHostToDevice);
        cudaMemcpy(ucand_dev,ucand.getPixels(), sizeof(float)*bsize*bsize, cudaMemcpyHostToDevice);
        ComputeCosts_noFeature_kernel<<<blocksPerGrid, threadsPerBlock>>>(ucand_dev, src.width(), src.height(), thrust::raw_pointer_cast(&candidates_dev[0]) , nsize*(rs+DMIR), dem_vars_dev, usr_var_dev,bsize);
        cudaThreadSynchronize();
        e_tmp = clock();
        match_time+=mstimer(s_tmp,e_tmp);

        s_tmp = clock();
        thrust::sort(candidates_dev.begin(),candidates_dev.end(),comp);
        cudaMemcpy(&candidates[0], thrust::raw_pointer_cast(&candidates_dev[0]), sizeof(cost_t)*nsize*(rs+DMIR), cudaMemcpyDeviceToHost);
        e_tmp = clock();
        sort_time+=mstimer(s_tmp,e_tmp);
        //cout<<"hi1\n";

        s_tmp = clock();
        Image patch=findCand(dest,candidates,src.width(), src.height(), prev,bsize,dx,dy);
        e_tmp = clock();
        find_time+=mstimer(s_tmp,e_tmp);
        //cout<<"hi2\n";

        s_tmp = clock();
        patch_merging(&dest, &patch, dx, dy,1,bsize/10.);
        e_tmp = clock();
        paste_time+=mstimer(s_tmp,e_tmp);

        //cout<<"hi3\n";

        cnum++;

        for (int i=0; i<bsize ; i++)	for (int j=0; j<bsize ; j++)
                if (xmaxi+i>=0 && ymaxi+j>=0 && xmaxi+i<dest.width() && ymaxi+j<dest.height() && dest(xmaxi+i,ymaxi+j)<=BG)
                    conf[(xmaxi+i)+(ymaxi+j)*dest.width()] = conf[xmaxi+ymaxi*dest.width()];
    }

    delete [] conf;

    bool havdata = false;
    bool empdata = false;
    for (int id = 0; id<dest.height()*dest.width(); id++)
        if (dest.getPixels()[id]>BG)
        {
            havdata = true;
            break;
        }
        else empdata = true;


    if (cnum==0 || empdata)
    {

        //cout<<"Number of targets: "<<cnum<<endl;
        while (true)
        {
            bool finish = true;
            for (int y = 0; y<dest.height(); y++)
                for (int x = 0; x<dest.width(); x++)
                    if (  ((!havdata) && x==0 && y==0 && dest(x,y)<=BG) ||  (dest(x,y)>BG && onBoundary(dest.getPixels(),dest.width(),dest.height(),x,y)))
                    {
                        finish = false;
                        dx = x-bsize/2;
                        dy = y-bsize/2;

                        s_tmp = clock();

                        vector<float> usr_var = noise_variances(tar_pyr,dx,dy,bsize);
                        Image ucand = dest.get_crop(dx,dy,dx+bsize-1,dy+bsize-1);

                        s_tmp = clock();
                        cudaMemcpy(usr_var_dev, &usr_var[0], sizeof(float)*NLEVEL, cudaMemcpyHostToDevice);
                        cudaMemcpy(ucand_dev,ucand.getPixels(), sizeof(float)*bsize*bsize, cudaMemcpyHostToDevice);
                        ComputeCosts_noFeature_kernel<<<blocksPerGrid, threadsPerBlock>>>(ucand_dev, src.width(), src.height(), thrust::raw_pointer_cast(&candidates_dev[0]) , nsize*(rs+DMIR), dem_vars_dev, usr_var_dev, bsize);
                        cudaThreadSynchronize();
                        e_tmp = clock();
                        match_time+=mstimer(s_tmp,e_tmp);

                        s_tmp = clock();
                        thrust::sort(candidates_dev.begin(),candidates_dev.end(),comp);
                        cudaMemcpy(&candidates[0], thrust::raw_pointer_cast(&candidates_dev[0]), sizeof(cost_t)*nsize*(rs+DMIR), cudaMemcpyDeviceToHost);
                        e_tmp = clock();
                        sort_time+=mstimer(s_tmp,e_tmp);

                        s_tmp = clock();
                        Image patch=findCand(dest,candidates,src.width(), src.height(), prev,bsize,dx,dy);
                        e_tmp = clock();
                        find_time+=mstimer(s_tmp,e_tmp);

                        s_tmp = clock();
                        patch_merging(&dest, &patch, dx, dy,1,bsize/10.);
                        e_tmp = clock();
                        paste_time+=mstimer(s_tmp,e_tmp);

                        cnum++;
                    }

            if (finish)	break;
        }

    }

    cudaUnbindTexture(utar_tex);
    cudaUnbindTexture(src_tex);

    delete [] dem_vars;

    cudaFree(usr_var_dev);
    cudaFree(ucand_dev);

    cudaFree(dnodes_dev);
    cudaFree(dem_vars_dev);

    cudaFreeArray(target_dev);
    cudaFreeArray(src_dev);

    cerr<<"\n\n*********** Non Feature matching GPU 2*******************\n";
    print_times();
    cerr<<" Candidates set: "<<cand_time<<"s\n";
    cerr<<" Number  of targets: "<<cnum<<"\n";
    cerr<<" Number  of cands: "<<nsize*(rs+DMIR)<<"\n";
    {
        end_t = clock();
        float elapsed = ((float)( end_t - start_t )) /CLOCKS_PER_SEC;
        cerr<<" Non-feature matching elapsed time: "<<elapsed<<" s.\n";
    }
    cerr<<"*********** End *******************\n\n";

}

__device__ float getCost_Feature(float* dest, float* dem_vars, node_t* dem_leafs, float* usr_var, node_t* uleafs, cost_t dem, int sw, int sh, int bsize, int dx, int dy)
{
    float tmp = 0.;

    node_t p = dem.org;

    if (use_noisestat)  tmp+=0.001*compare_variances(usr_var,dem_vars,dem.vpos,NLEVEL);

    if (use_angle) tmp+= 2*getDiffAng(uleafs,dem_leafs,dem.lpos,bsize,dem.lsize);

    if (use_profile)  tmp+= 5*get_diff_profile(sw,sh,bsize,p,dem.rot/DROT,uleafs,dem.lsize,dx,dy);
    if (use_ssd) tmp+= ssdf(dest, sw,sh,bsize,p.x,p.y,dem.rot/DROT);

    return tmp;
}

__global__ void ComputeCosts_Feature_kernel(float* dest, int sw, int sh, cost_t* candidates, int csize, float* dem_vars, node_t* dem_leafs, float* usr_var, node_t* uleafs, int bsize, int lsize, int dx, int dy)
{
    int k = blockDim.x * blockIdx.x + threadIdx.x;
    if (k<csize)
    {
        if ( (!candidates[k].skip) && lsize==candidates[k].lsize)
        {
            candidates[k].cost = getCost_Feature(dest, dem_vars, dem_leafs, usr_var, uleafs, candidates[k], sw, sh, bsize, dx, dy);
        }
        else if (!candidates[k].skip){
            candidates[k].cost = (INF/2)+ssdf(dest,sw,sh,bsize,candidates[k].org.x,candidates[k].org.y,candidates[k].rot/DROT);
        }
        else
        {
            candidates[k].cost = INF;
        }
    }
}

__global__ void ComputeCosts_noFeature_kernel(float* dest,float* bigtex, int ncand, cost_t* candidates, int csize, float* dem_vars, float* usr_var, int bsize)
{
    unsigned int k = blockDim.x * blockIdx.x + threadIdx.x;
    /*__shared__ float susr_var[BLOCK_SIZE];

    if ()

    __syncthreads;*/

    if (k<csize)
    {

        if ( (!candidates[k].skip))
        {
            candidates[k].cost = getCost_noFeature(dest, bigtex, ncand, dem_vars, usr_var, candidates[k], bsize);
        }
        else
        {
            candidates[k].cost = INF;
        }
    }
}

__global__ void ComputeCosts_Feature_kernel(float* dest,float* target, float* bigtex, int twidth, int ncand, cost_t* candidates, int csize, float* dem_vars, node_t* dem_leafs, float* usr_var, node_t* uleafs, int bsize, int lsize, int dx, int dy)
{
    unsigned int k = blockDim.x * blockIdx.x + threadIdx.x;

    if (k<csize)
    {

        if ( (!candidates[k].skip) && lsize==candidates[k].lsize)
        {
            candidates[k].cost = getCost_Feature(dest, target, bigtex, twidth, ncand, dem_vars, dem_leafs, usr_var, uleafs, candidates[k], bsize, dx,dy);
        }
        else if (!candidates[k].skip){
            candidates[k].cost = (INF/2)+ssdf(dest,bigtex,ncand,bsize,candidates[k].tnode.x);
        }
        else
        {
            candidates[k].cost = INF;
        }
    }
}

__global__ void buildCandidates_Feature_kernel(cost_t* candidates, node_t* dnodes,  node_t* dem_lsizes, int src_w, int src_h, int nsize,  int bsize)
{

    unsigned int k = blockDim.x * blockIdx.x + threadIdx.x;
    int rs = (360/DROT+DMIR);
    if (k<nsize*rs)
    {
        int kn = k%(nsize);
        int rt = k/(nsize);
        node_t cnode = dnodes[kn];

        int rx = cnode.x-bsize/2;
        int ry = cnode.y-bsize/2;
        int rot = rt*DROT;


        cost_t c(node_t(rx,ry),0,rot,0,0);
        {
            int kpos = kn+rt*nsize;

            c.lpos = dem_lsizes[kpos].x;
            c.lsize = dem_lsizes[kpos].y;
            c.vpos = kn*3;
            c.skip  =false;


            for (int id=0; id<bsize*bsize; id++)
            {


                int i = id%bsize;
                int j = id/bsize;
                float candv = getH_checkBounds(src_w, src_h, bsize, rx, ry, rt, i, j);
                if (candv<=BG)
                {
                    c.skip = true;
                }
            }

            candidates[kpos]=c;



        }
    }
}
__global__ void buildCandidates_Feature_kernel(cost_t* candidates, float* bigtex, float* src_ptr, node_t* dnodes, node_t* dem_lsizes, int nsize, int ncand, int src_w, int src_h, int bsize)
{
    unsigned int k = blockDim.x * blockIdx.x + threadIdx.x;
    int rs = (360/DROT+ DMIR);
    if (k<(nsize*rs))
    {

        //cout<<"Prepare candidates!\n";
        int cx=0, cy=0;

        int kn = k%(nsize);
        int rt = k/(nsize);

        node_t pnode = dnodes[kn];
        int rx = pnode.x-bsize/2;
        int ry = pnode.y-bsize/2;

        int rot = rt*DROT;



        cost_t c(node_t(rx,ry),0,rot,0,0);

        {
            c.vpos = kn*3;
            c.lpos = dem_lsizes[k].x;
            c.lsize = dem_lsizes[k].y;

            c.skip = false;

            cx = k;
            //cx = k*bsize*bsize;
            c.tnode = node_t(cx,cy);

            for (int id=0; id<bsize*bsize; id++)
            {

                int i = id%bsize;
                int j = id/bsize;
                //cout<<i<<"/"<j<<": "<<<<endl;
                float candv=getH(src_ptr,src_w,src_h,bsize,rx,ry,rt,i,j);
                (bigtex)[(cx)+(i+j*bsize)*ncand] = candv;

                if (candv<=BG)
                {
                    c.skip = true;
                }
            }
            candidates[k]=c;
        }
    }
}

__global__ void buildCandidates_noFeature_kernel(cost_t* candidates, float* bigtex, float* src_ptr, node_t* dnodes, int nsize, int ncand, int src_w, int src_h, int bsize)
{
    unsigned int k = blockDim.x * blockIdx.x + threadIdx.x;
    int rs = (360/DROT+ DMIR);
    if (k<(nsize*rs))
    {

        //cout<<"Prepare candidates!\n";
        int cx=0, cy=0;

        int kn = k%(nsize);
        int rt = k/(nsize);

        node_t pnode = dnodes[kn];
        int rx = pnode.x-bsize/2;
        int ry = pnode.y-bsize/2;

        int rot = rt*DROT;



        cost_t c(node_t(rx,ry),0,rot,0,0);

        {
            c.vpos = kn*3;

            c.skip = false;

            cx = k;
            //cx = k*bsize*bsize;
            c.tnode = node_t(cx,cy);

            for (int id=0; id<bsize*bsize; id++)
            {

                int i = id%bsize;
                int j = id/bsize;
                //cout<<i<<"/"<j<<": "<<<<endl;
                float candv=getH(src_ptr,src_w,src_h,bsize,rx,ry,rt,i,j);
                (bigtex)[(cx)+(i+j*bsize)*ncand] = candv;

                if (candv<=BG)
                {
                    c.skip = true;
                }
            }
            candidates[k]=c;
        }
    }
}

void match_Feature_gpu1(Terrain& dest, Tree& usr_features, Tree& dem_features, vector<Image>& tar_pyr, vector<Image>&  src_pyr,  int bsize)
{
    clock_t start_t, end_t, s_tmp, e_tmp;
    start_t = clock();

    match_time = 0;
    paste_time = 0;
    get_target = 0;
    find_time = 0;
    sort_time = 0;

    Image target = usr_features.msource;
    Image src = dem_features.msource;
    node_list dem_nodes = dem_features.processNodes;


    int nsize = dem_nodes.size();
    int rs = (360/DROT);

    int bigwidth  = nsize*bsize*(rs+DMIR)*bsize;
    int bigheight = 1;
    Image big(bigwidth,bigheight);
    float* bigtex = big.getPixels();
    //float* src_ptr = src.getPixels();
    int ncand = nsize*(rs+DMIR);
    //int ncand = 1;

    vector<cost_t> candidates(nsize*(rs+DMIR));
    node_t* dnodes = new node_t [nsize];

    s_tmp = clock();

    node_t* dem_leafs = new node_t [nsize*(rs+DMIR)*10];
    node_t* dem_lsizes = new node_t [nsize*(rs+DMIR)];

    s_tmp = clock();

    float* dem_vars =  matchingPrepocessing(dem_features, src_pyr, dem_nodes, dnodes, dem_leafs, dem_lsizes, bsize);
    {
        e_tmp = clock();
        float elapsed = ((float)( e_tmp - s_tmp )) /CLOCKS_PER_SEC;
        cerr<<"Get noise stats and control points elapsed time: "<<elapsed<<" s.\n";
        s_tmp = clock();
    }


    float* bigtex_dev;
    cudaMalloc((void**) &bigtex_dev,sizeof(float)*bigwidth*bigheight);
    float* src_dev;
    cudaMalloc((void**)&src_dev,sizeof(float)*src.width()*src.height());
    node_t* dnodes_dev;
    cudaMalloc((void**) &dnodes_dev, sizeof(node_t)*nsize);
    float* dem_vars_dev;
    cudaMalloc((void**) &dem_vars_dev, sizeof(float)*(nsize*NLEVEL));


    cudaMemcpy(src_dev, src.getPixels(), sizeof(float)*src.width()*src.height(), cudaMemcpyHostToDevice);
    cudaMemcpy(dnodes_dev, dnodes, sizeof(node_t)*nsize, cudaMemcpyHostToDevice);
    cudaMemcpy(dem_vars_dev, dem_vars, sizeof(float)*(nsize*NLEVEL), cudaMemcpyHostToDevice);


    thrust::device_vector<cost_t> candidates_dev (nsize*(rs+DMIR));

    node_t* dem_lsizes_dev;
    cudaMalloc((void**) &dem_lsizes_dev, sizeof(node_t)*nsize*(rs+DMIR));
    node_t* dem_leafs_dev;
    cudaMalloc((void**) &dem_leafs_dev, sizeof(node_t)*nsize*(rs+DMIR)*10);

    cudaMemcpy(dem_lsizes_dev, dem_lsizes, sizeof(node_t)*nsize*(rs+DMIR), cudaMemcpyHostToDevice);
    cudaMemcpy(dem_leafs_dev, dem_leafs, sizeof(node_t)*nsize*(rs+DMIR)*10, cudaMemcpyHostToDevice);

    float* target_dev;
    cudaMalloc((void**) &target_dev, sizeof(float)*target.width()*target.height());
    cudaMemcpy(target_dev, target.getPixels(), sizeof(float)*target.width()*target.height(), cudaMemcpyHostToDevice);

    float cand_time=0;
    int threadsPerBlock = BLOCK_SIZE;
    int blocksPerGrid =  ((nsize*(rs+DMIR)) / threadsPerBlock)+1;

    //buildCandidates_Feature(&candidates[0], bigtex, src_ptr, dnodes, dem_lsizes,nsize, ncand, src.getWidth(),src.getHeight(),bsize);
    buildCandidates_Feature_kernel<<<blocksPerGrid, threadsPerBlock>>>(thrust::raw_pointer_cast(&candidates_dev[0]), bigtex_dev, src_dev, dnodes_dev, dem_lsizes_dev,nsize, ncand, src.getWidth(),src.getHeight(),bsize);
    cudaMemcpy(bigtex, bigtex_dev, sizeof(float)*bigwidth*bigheight, cudaMemcpyDeviceToHost);
    //cudaMemcpy(bigtex_dev, bigtex, sizeof(float)*nsize*bsize*(rs+DMIR)*bsize, cudaMemcpyHostToDevice);
    //cudaMemcpy(thrust::raw_pointer_cast(&candidates_dev[0]),&candidates[0], sizeof(cost_t)*nsize*(rs+DMIR), cudaMemcpyHostToDevice);

    {
        e_tmp = clock();
        float elapsed = ((float)( e_tmp - s_tmp )) /CLOCKS_PER_SEC;
        cand_time += elapsed;
    }

    //big.savePGM("bigtmp_Feature.pgm");
    //cudaMemcpy(src_ptr, src_dev, sizeof(float)*src.width()*src.height(), cudaMemcpyDeviceToHost);
    //Image cand(bsize, bsize);
    //getCand(cand, bigtex, ncand, 4*bsize*bsize);

    //cand.savePGM("cand.pgm");


    node_list usr_nodes = usr_features.processNodes;
    cost_t prev (node_t(-1,-1),0,0,0,0);

    float* usr_var_dev;
    cudaMalloc((void**) &usr_var_dev,sizeof(float)*NLEVEL);
    node_t* uleafs_dev;
    cudaMalloc((void**) &uleafs_dev,sizeof(node_t)*10);
    float* ucand_dev;
    cudaMalloc(&ucand_dev, sizeof(float)*bsize*bsize);

    cout<<"Start matching feature patches! "<<usr_nodes.size()<<" from "<<nsize*(rs+DMIR)<<" candidates\n";



    int dx,dy;
    for (node_list::const_iterator it = usr_nodes.begin(); it != usr_nodes.end(); it++ )
    {

        s_tmp = clock();
        node_t cnode = *it;
        dx = cnode.x-bsize/2;
        dy = cnode.y-bsize/2;

        vector<node_t> uleafs = getChildren(cnode,usr_features,bsize);
        vector<float> usr_var =  noise_variances(tar_pyr,dx,dy,bsize);
        Image ucand = dest.get_crop(dx,dy,dx+bsize-1,dy+bsize-1);
        int lsize = min(uleafs.size(),10);

        s_tmp = clock();
        cudaMemcpy(ucand_dev, ucand.getPixels(), sizeof(float)*bsize*bsize, cudaMemcpyHostToDevice);
        cudaMemcpy(uleafs_dev, &uleafs[0], sizeof(node_t)*lsize, cudaMemcpyHostToDevice);
        cudaMemcpy(usr_var_dev, &usr_var[0], sizeof(float)*NLEVEL, cudaMemcpyHostToDevice);
        ComputeCosts_Feature_kernel<<<blocksPerGrid, threadsPerBlock>>>(ucand_dev,target_dev,bigtex_dev,target.width(),ncand,thrust::raw_pointer_cast(&candidates_dev[0]),nsize*(rs+DMIR),dem_vars_dev,dem_leafs_dev,usr_var_dev,uleafs_dev,bsize, uleafs.size(), dx,dy);
        cudaThreadSynchronize();
        e_tmp = clock();
        match_time+=mstimer(s_tmp,e_tmp);

        s_tmp = clock();
        thrust::sort(candidates_dev.begin(),candidates_dev.end(),comp);
        cudaMemcpy(&candidates[0], thrust::raw_pointer_cast(&candidates_dev[0]), sizeof(cost_t)*nsize*(rs+DMIR), cudaMemcpyDeviceToHost);
        e_tmp = clock();
        sort_time+=mstimer(s_tmp,e_tmp);

        s_tmp = clock();
        Image patch = findCand(dest, bigtex, ncand, candidates,prev, bsize, dx,dy);
        e_tmp = clock();
        find_time+=mstimer(s_tmp,e_tmp);

        s_tmp = clock();
        patch_merging(&dest, &patch, dx, dy,1,bsize/10.);
        e_tmp = clock();
        paste_time+=mstimer(s_tmp,e_tmp);

        //dest.savePGM("tmp.pgm");         cin.get();

    }

    delete [] dem_vars;
    delete [] dem_leafs;
    delete [] dnodes;
    delete [] dem_lsizes;
    cudaFree(dem_lsizes_dev);

    cudaFree(usr_var_dev);
    cudaFree(uleafs_dev);
    cudaFree(ucand_dev);

    cudaFree(dem_vars_dev);
    cudaFree(dem_leafs_dev);
    cudaFree(dnodes_dev);

    cudaFree(src_dev);
    cudaFree(target_dev);
    cudaFree(bigtex_dev);

    cerr<<"\n\n*********** Feature matching GPU 1*******************\n";
    print_times();
    cerr<<" Candidates set: "<<cand_time<<"s\n";
    cerr<<" Number  of targets: "<<usr_nodes.size()<<"\n";
    cerr<<" Number  of cands: "<<nsize*(rs+DMIR)<<"\n";
    {
        end_t = clock();
        float elapsed = ((float)( end_t - start_t )) /CLOCKS_PER_SEC;
        cerr<<" Feature matching elapsed time: "<<elapsed<<" s.\n";
    }
    cerr<<"*********** End *******************\n\n";
}


void match_Feature_gpu2(Terrain& dest, Tree& usr_features, Tree& dem_features, vector<Image>& tar_pyr, vector<Image>&  src_pyr, int bsize)
{
    clock_t start_t, end_t, s_tmp, e_tmp;
    start_t = clock();

    match_time = 0;
    paste_time = 0;
    get_target = 0;
    find_time = 0;
    sort_time = 0;

    Image& target = usr_features.msource;
    Image& src = dem_features.msource;
    node_list dem_nodes = dem_features.processNodes;



    int nsize = dem_nodes.size();
    int rs = 360/DROT;

    node_t* dem_leafs = new node_t [nsize*(rs+DMIR)*10];
    node_t* dem_lsizes = new node_t [nsize*(rs+DMIR)];
    node_t* dnodes = new node_t [nsize];


    s_tmp = clock();

    float* dem_vars =  matchingPrepocessing(dem_features, src_pyr, dem_nodes, dnodes, dem_leafs, dem_lsizes, bsize);

    {
        e_tmp = clock();
        float elapsed = ((float)( e_tmp - s_tmp )) /CLOCKS_PER_SEC;
        cerr<<"Get noise stats and control points elapsed time: "<<elapsed<<" s.\n";
        s_tmp = clock();
    }

    float cand_time = 0;
    //cout<<"Prepare candidates! "<< nsize<<"\n";
    s_tmp = clock();

    float* dem_vars_dev;
    cudaMalloc((void**) &dem_vars_dev, sizeof(float)*nsize*3);
    thrust::device_vector<cost_t> candidates_dev(nsize*(rs+DMIR));
    node_t* dnodes_dev;
    cudaMalloc((void**) &dnodes_dev, sizeof(node_t)*nsize);
    node_t* dem_lsizes_dev;
    cudaMalloc((void**) &dem_lsizes_dev, sizeof(node_t)*nsize*(rs+DMIR));
    node_t* dem_leafs_dev;
    cudaMalloc((void**) &dem_leafs_dev, sizeof(node_t)*nsize*(rs+DMIR)*10);

    cudaMemcpy(dem_vars_dev, dem_vars, sizeof(float)*nsize*3, cudaMemcpyHostToDevice);
    cudaMemcpy(dnodes_dev, dnodes, sizeof(node_t)*nsize, cudaMemcpyHostToDevice);
    cudaMemcpy(dem_lsizes_dev, dem_lsizes, sizeof(node_t)*nsize*(rs+DMIR), cudaMemcpyHostToDevice);
    cudaMemcpy(dem_leafs_dev, dem_leafs, sizeof(node_t)*nsize*(rs+DMIR)*10, cudaMemcpyHostToDevice);



    int threadsPerBlock = BLOCK_SIZE;
    int blocksPerGrid =  ((nsize*(rs+DMIR)) / threadsPerBlock)+1;

    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
    cudaArray* src_dev;
    cudaMallocArray(&src_dev, &channelDesc, src.width(), src.height());
    cudaMemcpyToArray(src_dev, 0, 0, src.getPixels(), sizeof(float)*src.width()*src.height(), cudaMemcpyHostToDevice);
    // Set texture parameters
    src_tex.addressMode[0] = cudaAddressModeClamp;
    src_tex.addressMode[1] = cudaAddressModeClamp;
    src_tex.filterMode = cudaFilterModePoint;
    src_tex.normalized = false;
    cudaBindTextureToArray(src_tex, src_dev, channelDesc);

    cudaArray* target_dev;
    cudaMallocArray(&target_dev, &channelDesc, target.width(), target.height());
    cudaMemcpyToArray(target_dev, 0, 0, target.getPixels(), sizeof(float)*target.width()*target.height(), cudaMemcpyHostToDevice);
    // Set texture parameters
    utar_tex.addressMode[0] = cudaAddressModeClamp;
    utar_tex.addressMode[1] = cudaAddressModeClamp;
    utar_tex.filterMode = cudaFilterModePoint;
    utar_tex.normalized = false;
    cudaBindTextureToArray(utar_tex, target_dev, channelDesc);

    {
        buildCandidates_Feature_kernel<<<blocksPerGrid, threadsPerBlock>>>(thrust::raw_pointer_cast(&candidates_dev[0]), dnodes_dev, dem_lsizes_dev,  src.width(), src.height(), nsize, bsize);
        cudaThreadSynchronize();
        e_tmp = clock();
        float elapsed = ((float)( e_tmp - s_tmp )) /CLOCKS_PER_SEC;
        cand_time += elapsed;
    }

    delete [] dem_vars;
    delete [] dem_leafs;
    delete [] dnodes;
    delete [] dem_lsizes;
    cudaFree(dem_lsizes_dev);

    node_list usr_nodes = usr_features.processNodes;
    cost_t prev(node_t(-1,-1),0,0,0,0);

    float* usr_var_dev;
    cudaMalloc((void**) &usr_var_dev,sizeof(float)*NLEVEL);
    node_t* uleafs_dev;
    cudaMalloc((void**) &uleafs_dev,sizeof(node_t)*10);

    float* ucand_dev;
    cudaMalloc(&ucand_dev, sizeof(float)*bsize*bsize);

    vector<cost_t> candidates (nsize*(rs+DMIR));


    cout<<"Start matching feature patches! "<<usr_nodes.size()<<" from "<<nsize*(rs+DMIR)<<" candidates\n";

    int dx,dy;
    for (node_list::const_iterator it = usr_nodes.begin(); it != usr_nodes.end(); it++ )
    {
        s_tmp = clock();

        node_t cnode = *it;
        dx = cnode.x-bsize/2;
        dy = cnode.y-bsize/2;

        vector<node_t> uleafs = getChildren(cnode,usr_features,bsize);
        vector<float> usr_var =  noise_variances(tar_pyr,dx,dy,bsize);
        int lsize = min(uleafs.size(),10);

        Image ucand = dest.get_crop(dx,dy,dx+bsize-1,dy+bsize-1);

        s_tmp = clock();
        cudaMemcpy(ucand_dev, ucand.getPixels(), sizeof(float)*bsize*bsize, cudaMemcpyHostToDevice);
        cudaMemcpy(uleafs_dev, &uleafs[0], sizeof(node_t)*lsize, cudaMemcpyHostToDevice);
        cudaMemcpy(usr_var_dev, &usr_var[0], sizeof(float)*NLEVEL, cudaMemcpyHostToDevice);
        ComputeCosts_Feature_kernel<<<blocksPerGrid, threadsPerBlock>>>(ucand_dev, src.width(), src.height(), thrust::raw_pointer_cast(&candidates_dev[0]) , nsize*(rs+DMIR), dem_vars_dev, dem_leafs_dev, usr_var_dev, uleafs_dev,bsize, lsize, dx, dy);
        cudaThreadSynchronize();
        e_tmp = clock();
        match_time+=mstimer(s_tmp,e_tmp);

        s_tmp = clock();
        thrust::sort(candidates_dev.begin(),candidates_dev.end(),comp);
        cudaMemcpy(&candidates[0], thrust::raw_pointer_cast(&candidates_dev[0]), sizeof(cost_t)*nsize*(rs+DMIR), cudaMemcpyDeviceToHost);
        e_tmp = clock();
        sort_time+=mstimer(s_tmp,e_tmp);

        s_tmp = clock();
        Image patch=findCand(dest,candidates, src.width(), src.height(), prev,bsize,dx,dy);
        e_tmp = clock();
        find_time+=mstimer(s_tmp,e_tmp);

        s_tmp = clock();
        patch_merging(&dest, &patch, dx, dy,1,bsize/10.);
        e_tmp = clock();
        paste_time+=mstimer(s_tmp,e_tmp);
    }

    cudaUnbindTexture(src_tex);
    cudaUnbindTexture(utar_tex);

    cudaFree(usr_var_dev);
    cudaFree(uleafs_dev);
    cudaFree(ucand_dev);
    //cudaFreeArray(ucand_dev);

    cudaFree(dnodes_dev);

    //cudaFree(candidates_dev);
    //delete prev;
    cudaFree(dem_vars_dev);
    cudaFree(dem_leafs_dev);
    cudaFree(dnodes_dev);
    //cudaFree(src_dev);
    cudaFreeArray(src_dev);
    cudaFreeArray(target_dev);
    cerr<<"\n\n*********** Feature matching GPU 2*******************\n";
    print_times();
    cerr<<" Candidates set: "<<cand_time<<"s\n";
    cerr<<" Number  of targets: "<<usr_nodes.size()<<"\n";
    cerr<<" Number  of cands: "<<nsize*(rs+DMIR)<<"\n";
    {
        end_t = clock();
        float elapsed = ((float)( end_t - start_t )) /CLOCKS_PER_SEC;
        cerr<<" Feature matching elapsed time: "<<elapsed<<" s.\n";
    }
    cerr<<"*********** End *******************\n\n";
}

void match_noFeature_gpu1(Terrain& dest, Image& src, Image& target, node_list dem_nodes, vector<Image>& tar_pyr, vector<Image>&  src_pyr,  int bsize, int osize)
{
    clock_t start_t, end_t, s_tmp, e_tmp;
    start_t = clock();

    match_time = 0;
    paste_time = 0;
    get_target = 0;
    find_time = 0;
    sort_time = 0;

    int nsize = dem_nodes.size();
    int rs = (360/DROT);

    int bigwidth  = nsize*bsize*(rs+DMIR)*bsize;
    int bigheight = 1;
    Image big(bigwidth,bigheight);
    float* bigtex = big.getPixels();
    //float* src_ptr = src.getPixels();
    int ncand = nsize*(rs+DMIR);
    //int ncand = 1;

    vector<cost_t> candidates(nsize*(rs+DMIR));
    node_t* dnodes = new node_t [nsize];

    s_tmp = clock();

    float* dem_vars =  matchingPrepocessing(src_pyr, dem_nodes, dnodes, bsize);

    {
        e_tmp = clock();
        float elapsed = ((float)( e_tmp - s_tmp )) /CLOCKS_PER_SEC;
        cerr<<"Get noise stats and control points elapsed time: "<<elapsed<<" s.\n";
        s_tmp = clock();
    }


    float* bigtex_dev;
    cudaMalloc((void**) &bigtex_dev,sizeof(float)*bigwidth*bigheight);
    float* src_dev;
    cudaMalloc((void**)&src_dev,sizeof(float)*src.width()*src.height());
    node_t* dnodes_dev;
    cudaMalloc((void**) &dnodes_dev, sizeof(node_t)*nsize);
    float* dem_vars_dev;
    cudaMalloc((void**) &dem_vars_dev, sizeof(float)*(nsize*NLEVEL));


    cudaMemcpy(src_dev, src.getPixels(), sizeof(float)*src.width()*src.height(), cudaMemcpyHostToDevice);
    cudaMemcpy(dnodes_dev, dnodes, sizeof(node_t)*nsize, cudaMemcpyHostToDevice);
    cudaMemcpy(dem_vars_dev, dem_vars, sizeof(float)*(nsize*NLEVEL), cudaMemcpyHostToDevice);


    thrust::device_vector<cost_t> candidates_dev (nsize*(rs+DMIR));

    float cand_time=0;
    int threadsPerBlock = BLOCK_SIZE;
    int blocksPerGrid =  ((nsize*(rs+DMIR)) / threadsPerBlock)+1;
    buildCandidates_noFeature_kernel<<<blocksPerGrid, threadsPerBlock>>>(thrust::raw_pointer_cast(&candidates_dev[0]), bigtex_dev, src_dev, dnodes_dev, nsize, ncand, src.getWidth(),src.getHeight(),bsize);
    cudaMemcpy(bigtex, bigtex_dev, sizeof(float)*bigwidth*bigheight, cudaMemcpyDeviceToHost);

    {
        e_tmp = clock();
        float elapsed = ((float)( e_tmp - s_tmp )) /CLOCKS_PER_SEC;
        cand_time += elapsed;
    }

    cost_t prev (node_t(-1,-1),0,0,0,0);

    int cnum = nsize*(rs+DMIR);
    cout<<"Start matching non-feature patches! from "<<cnum<<"\n";

    float* usr_var_dev;
    cudaMalloc((void**) &usr_var_dev,sizeof(float)*NLEVEL);
    float* ucand_dev;
    cudaMalloc((void **) &ucand_dev, sizeof(float)*bsize*bsize);
    cnum = 0;

    int dx = -5*bsize, dy=-5*bsize;

    list<cost_t> omega;
    float* conf = new float [dest.width()*dest.height()];
    for (int x=0; x<dest.width(); x++)	for (int y=0; y<dest.height(); y++)
            if (dest(x,y)>BG)
                conf[x+y*dest.width()] = 1;
            else
                conf[x+y*dest.width()] = 0;
    while (true)
    {
        s_tmp = clock();

        getNextTarget (dest.getPixels(), conf, omega,dest.width(), dest.height(), bsize, dx+bsize/2, dy+bsize/2);
        if (omega.size()==0) break;


        int xmaxi = -1;
        int ymaxi = -1;
        float maxv = -1e6f;
        float tmp=0.0f;

        for (list<cost_t>::iterator it=omega.begin(); it!=omega.end(); it++)
        {
            tmp = (*it).cost;
            if (tmp>maxv+1e-6f)
            {
                node_t p = (*it).org;
                maxv = tmp;
                xmaxi = p.x;
                ymaxi = p.y;
            }
        }

        e_tmp = clock();
        get_target+=mstimer(s_tmp,e_tmp);

        dx=(xmaxi)-bsize/2;
        dy=(ymaxi)-bsize/2;


        vector<float> usr_var = noise_variances(tar_pyr,dx,dy,bsize);
        Image ucand = dest.get_crop(dx,dy,dx+bsize-1,dy+bsize-1);

        s_tmp = clock();
        cudaMemcpy(usr_var_dev, &usr_var[0], sizeof(float)*NLEVEL, cudaMemcpyHostToDevice);
        cudaMemcpy(ucand_dev,ucand.getPixels(), sizeof(float)*bsize*bsize, cudaMemcpyHostToDevice);
        ComputeCosts_noFeature_kernel<<<blocksPerGrid, threadsPerBlock>>>(ucand_dev, bigtex_dev, ncand, thrust::raw_pointer_cast(&candidates_dev[0]), nsize*(rs+DMIR), dem_vars_dev, usr_var_dev,bsize);
        cudaThreadSynchronize();
        e_tmp = clock();
        match_time+=mstimer(s_tmp,e_tmp);

        s_tmp = clock();
        thrust::sort(candidates_dev.begin(),candidates_dev.end(),comp);
        cudaMemcpy(&candidates[0], thrust::raw_pointer_cast(&candidates_dev[0]), sizeof(cost_t)*nsize*(rs+DMIR), cudaMemcpyDeviceToHost);
        e_tmp = clock();
        sort_time+=mstimer(s_tmp,e_tmp);

        s_tmp = clock();
        Image patch = findCand(dest, bigtex, ncand, candidates,prev, bsize, dx,dy);
        e_tmp = clock();
        find_time+=mstimer(s_tmp,e_tmp);

        s_tmp = clock();
        patch_merging(&dest, &patch, dx, dy,1,bsize/10.);
        e_tmp = clock();
        paste_time+=mstimer(s_tmp,e_tmp);

        cnum++;

        for (int i=0; i<bsize ; i++)	for (int j=0; j<bsize ; j++)
                if (xmaxi+i>=0 && ymaxi+j>=0 && xmaxi+i<dest.width() && ymaxi+j<dest.height() && dest(xmaxi+i,ymaxi+j)<=BG)
                    conf[(xmaxi+i)+(ymaxi+j)*dest.width()] = conf[xmaxi+ymaxi*dest.width()];
    }

    delete [] conf;

    bool havdata = false;
    bool empdata = false;
    for (int id = 0; id<dest.height()*dest.width(); id++)
        if (dest.getPixels()[id]>BG)
        {
            havdata = true;
            break;
        }
        else
            empdata = true;


    if (cnum==0 || empdata)
    {

        //cout<<"Number of targets: "<<cnum<<endl;
        while (true)
        {
            bool finish = true;
            for (int y = 0; y<dest.height(); y++)
                for (int x = 0; x<dest.width(); x++)
                    if (  ((!havdata) && x==0 && y==0 && dest(x,y)<=BG) ||  (dest(x,y)>BG && onBoundary(dest.getPixels(),dest.width(),dest.height(),x,y)))
                    {
                        finish = false;
                        dx = x-bsize/2;
                        dy = y-bsize/2;

                        s_tmp = clock();

                        vector<float> usr_var = noise_variances(tar_pyr,dx,dy,bsize);
                        Image ucand = dest.get_crop(dx,dy,dx+bsize-1,dy+bsize-1);

                        s_tmp = clock();
                        cudaMemcpy(usr_var_dev, &usr_var[0], sizeof(float)*NLEVEL, cudaMemcpyHostToDevice);
                        cudaMemcpy(ucand_dev,ucand.getPixels(), sizeof(float)*bsize*bsize, cudaMemcpyHostToDevice);
                        ComputeCosts_noFeature_kernel<<<blocksPerGrid, threadsPerBlock>>>(ucand_dev, bigtex_dev, ncand, thrust::raw_pointer_cast(&candidates_dev[0]), nsize*(rs+DMIR), dem_vars_dev, usr_var_dev,bsize);
                        cudaThreadSynchronize();
                        e_tmp = clock();
                        match_time+=mstimer(s_tmp,e_tmp);

                        s_tmp = clock();
                        thrust::sort(candidates_dev.begin(),candidates_dev.end(),comp);
                        cudaMemcpy(&candidates[0], thrust::raw_pointer_cast(&candidates_dev[0]), sizeof(cost_t)*nsize*(rs+DMIR), cudaMemcpyDeviceToHost);
                        e_tmp = clock();
                        sort_time+=mstimer(s_tmp,e_tmp);

                        s_tmp = clock();
                        Image patch = findCand(dest, bigtex, ncand, candidates,prev, bsize, dx,dy);
                        e_tmp = clock();
                        find_time+=mstimer(s_tmp,e_tmp);

                        s_tmp = clock();
                        patch_merging(&dest, &patch, dx, dy,1,bsize/10.);
                        e_tmp = clock();
                        paste_time+=mstimer(s_tmp,e_tmp);

                        cnum++;
                    }

            if (finish)	break;
        }

    }

    cudaFree(usr_var_dev);
    cudaFree(ucand_dev);

    delete [] dnodes;
    delete [] dem_vars;

    cudaFree(bigtex_dev);
    cudaFree(src_dev);
    cudaFree(dnodes_dev);
    cudaFree(dem_vars_dev);

    cerr<<"\n\n*********** Non Feature matching GPU 1*******************\n";
    print_times();
    cerr<<" Candidates set: "<<cand_time<<"s\n";
    cerr<<" Number  of targets: "<<cnum<<"\n";
    cerr<<" Number  of cands: "<<nsize*(rs+DMIR)<<"\n";
    {
        end_t = clock();
        float elapsed = ((float)( end_t - start_t )) /CLOCKS_PER_SEC;
        cerr<<" Non-feature matching elapsed time: "<<elapsed<<" s.\n";
    }
    cerr<<"*********** End *******************\n\n";
}

void match_noFeature(Terrain& dest, Image& src, Image& target, node_list dem_nodes, vector<Image>& tar_pyr, vector<Image>&  src_pyr, int bsize, int osize)
{
    match_noFeature_gpu2(dest, src, target, dem_nodes, tar_pyr, src_pyr, bsize, osize);
}

void match_Feature(Terrain& dest, Tree& usr_features, Tree& dem_features, vector<Image>& tar_pyr, vector<Image>&  src_pyr, int bsize)
{
    match_Feature_gpu2(dest, usr_features, dem_features, tar_pyr, src_pyr, bsize);
}
