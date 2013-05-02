#pragma once

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <cfloat>

#define PI acosf(-1)
#define MAX_FLOAT FLT_MAX

#define vPAD 255.

#define vSINK 195.
#define vsSINK 127.
#define vsSOURCE 63
#define vSOURCE 31.

#define NLEVEL 3

#define DEFAULT_DIMX 5000
#define DEFAULT_DIMY 5000

#define MAX_VAL  255


#define DEGREES_PER_PIXEL 0.6f
#define RADIANS_PER_PIXEL 0.002f
#define UNITS_PER_PIXEL 10.2f
#define UNITS_PER_WHEELTICK 0.35f
#define ZOOM_FACTOR .04f
#define WATER_UNITS 5.6f

#define INF MAX_FLOAT
#define DEFAULT_DIMX 5000
#define DEFAULT_DIMY 5000

#define BG 0.0 //background color

// Includes
#include <iostream>

#include <map>
#include <list>
#include <vector>
#include <queue>
#include <algorithm>
#include <utility>
#include <fstream>
#include <sys/stat.h>

using namespace std;

//typedef unsigned int uint;
typedef pair<int,int> pair_t;   //t: public pair_t

class node_t {
    public:
        int x; int y;
        node_t(){x=0; y=0;}
        node_t(int _x, int _y) {x=_x; y=_y;}
        node_t(const node_t& v){
            x = v.x; y= v.y;
        }

        void operator ()(const int x0, const int y0) {	x= x0; y= y0;	}

        //! set to value
        const node_t &operator =(const node_t &v){
            x= v.x; y= v.y;
            return *this;
        }

        //! test for equality
        bool operator==(const node_t &v) const{	return (x==v.x && y==v.y);	}

        //! test for inequality
        bool operator!=(const node_t &v)  const{	return (x!=v.x || y!=v.y);	}

        bool operator<(const node_t &v)  const{	pair_t pair1(x,y); pair_t pair2(v.x,v.y); return (pair1<pair2);	}
};

class point_t{
    public:
        float x; float y;
        point_t(){x=0; y=0;}
        point_t(float _x) {x=_x; y=0;}
        point_t(float _x, float _y) {x=_x; y=_y;}

        point_t  & operator=(const point_t &rhs){
            if (this != &rhs) {
                x = rhs.x;
                y  =rhs.y;
            }
            return *this;
        }

         bool operator==(const point_t &other) const {
            return (x == other.x && y==other.y);
          }

        bool operator!=(const point_t &other) const {
            return !(*this == other);
          }
};

class cost_t{
    public:
        node_t org; int rot; int mir; int num; float cost;
        int lpos;  // leafs position
        int vpos;  //noise variance
        int lsize;
        bool skip;
        node_t cnode;
        node_t tnode;

        cost_t(){
        	org.x = 0;
        	org.y = 0;
        	cost = 0;
        }

        cost_t( node_t _org, float _cost=0., int _rot=0, int _mir=0, int _num=0) {
            org=_org; cost=_cost;  mir=_mir; rot = _rot; num=_num;
            lpos = 0;
            vpos = 0;
            lsize = 0;
            skip = false;
        }
        /*cost_t(int _k) {k=_k; cost=0; mir=0; rot = 0; num=0;}
        cost_t(int _k, int _cost) {k=_k; cost=_cost;  mir=0; rot = 0; num=0;}
        cost_t(int _k, int _cost, int _rot, int _mir, int _num=0) {k=_k; cost=_cost;  mir=_mir; rot = _rot; num=_num;}*/

        cost_t  & operator=(const cost_t &rhs){
            if (this != &rhs) {
                org = rhs.org;
                cost = rhs.cost;
                rot = rhs.rot;
                mir = rhs.mir;
                num = rhs.num;

                lpos = rhs.lpos;
                lsize = rhs.lsize;
                vpos = rhs.vpos;
                cnode = rhs.cnode;
                tnode = rhs.tnode;

                skip = rhs.skip;
            }
            return *this;
        }

        bool operator==(const cost_t &other) const {
            return (org == other.org && rot==other.rot);//  && mir==other.mir && rot==other.rot && num==other.num); //&& cost==other.cost
          }

         bool operator<(const cost_t &other) const {
            return (cost<other.cost);
          }

         bool operator!=(const cost_t &other) const {
            return !(*this == other);
          }
};



typedef std::vector<node_t> node_vec;
typedef std::list<node_t> node_list;
typedef std::map<node_t,node_list > map_t;

typedef std::pair<node_t, bool> status_t;
typedef std::vector< status_t > graphcut_t;
