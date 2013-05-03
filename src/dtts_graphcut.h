#include "external/maxflow/graph.h"
#include "dtts_image.h"
using namespace Dtts;

#define MAX_CAPACITY 999999999.

typedef Graph<float,float,float> GraphType;
typedef pair<node_t,bool> cut_node;


Image performCut(Image& dest, Image& psrc, int dx, int dy, float& cost){

	 vector<cut_node> cut;
     Image mask(psrc.width(),psrc.height());

     // Determine overlapping area
	for (int x=0; x<psrc.width(); x++) for (int y=0; y<psrc.height(); y++)
        if (x+dx<dest.width() && y+dy<dest.height() && dest(x+dx,y+dy)!=0.)
            cut.push_back( cut_node (node_t(x,y),true) );

     if (cut.size()==0){
         cost = 0.;	return mask;
     }

	 GraphType *G= new GraphType(cut.size(),4*cut.size());
	 G->add_node(cut.size());

	 for (unsigned int k=0; k<cut.size(); k++){
        node_t p = cut[k].first;

		for (unsigned int l=k+1; l<cut.size(); l++){
            node_t q = cut[l].first;
                if ( ( p.x==q.x || p.y==q.y ) && (p.x-q.x)*(p.x-q.x)+(p.y-q.y)*(p.y-q.y)==1 ){
                    float val = abs( psrc(p.x,p.y) - dest(p.x+dx,p.y+dy) )+abs( psrc(q.x,q.y) - dest(q.x+dx,q.y+dy) );
                    G->add_edge( k,l,val,val);
				}
            }


		bool bg=false;
		for(int m=-1; m<2; m++){
            if (bg)	break;
            for(int n=-1; n<2; n++){
                if(dest.atXY(p.x+dx+m,p.y+dy+n)==0.){
                    bg = true; break;
				}
            }
        }

        if (bg){ //Coming from source
            G->add_tweights(k,MAX_CAPACITY,0);
            continue;
        }

		bool fg=false;
        if ((p.y==0 || p.x==0 || p.y==psrc.height()-1 || p.x==psrc.width()-1)){
			int nx; int ny;

            nx = p.x; ny = p.y-1;
			if ( ny<0 && ny+dy>=0 && dest(dx+nx, dy+ny)!=0. )	fg = true;

			nx = p.x; ny = p.y+1;
			if ( ny>=psrc.height() && ny+dy<dest.height() && dest(dx+nx, dy+ny)!=0. )	fg = true;

			nx = p.x-1; ny = p.y;
			if ( nx<0 && nx+dx>=0 && dest(dx+nx, dy+ny)!=0. )	fg = true;

			nx = p.x+1; ny = p.y;
			if ( nx>=psrc.width() && nx+dx<dest.width() && dest(dx+nx, dy+ny)!=0. )	fg = true;

        }

        if (fg){ //Coming from sink
            G->add_tweights(k,0,MAX_CAPACITY);
            //mask(p.x,p.y) = vsSINK; //to comment
        }


	 }


	cost = G->maxflow();
	//cout<<"Flow: "<<f<<endl;

	for (unsigned int k=0; k<cut.size(); k++){
		node_t pnode = cut[k].first;
		if(G->what_segment(k) == GraphType::SINK){
			cut[k].second=false;
		}
	}
	delete G;

    //return mask;

    Image tmp(mask);
    if (cut.size()>0){

		for (unsigned int k=0; k<cut.size(); k++){
			 node_t p = cut[k].first;
			 if (!cut[k].second)   tmp(p.x,p.y) = vSINK;
			 else 	tmp(p.x,p.y) = vSOURCE;
        }

        for (int i=0; i<mask.width(); i++) for (int j=0; j<mask.height(); j++)
            if (tmp(i,j)==vSINK){
                mask(i,j)=vSINK;
                 if ( tmp.atXY(i-1,j)==vSOURCE || tmp.atXY(i+1,j)==vSOURCE || tmp.atXY(i,j-1)==vSOURCE || tmp.atXY(i,j+1)==vSOURCE ||
                      tmp.atXY(i-1,j-1)==vSOURCE || tmp.atXY(i+1,j-1)==vSOURCE || tmp.atXY(i-1,j+1)==vSOURCE ||tmp.atXY(i+1,j+1)==vSOURCE )
                        mask(i,j) = vsSINK;
            }
            else if (tmp(i,j)==vSOURCE){
                mask(i,j)=vSOURCE;
                if ( tmp.atXY(i-1,j)==vSINK || tmp.atXY(i+1,j)==vSINK || tmp.atXY(i,j-1)==vSINK || tmp.atXY(i,j+1)==vSINK ||
                      tmp.atXY(i-1,j-1)==vSINK || tmp.atXY(i+1,j-1)==vSINK || tmp.atXY(i-1,j+1)==vSINK ||tmp.atXY(i+1,j+1)==vSINK )
                        mask(i,j) = vsSOURCE;
            }

    }

	return mask;
}


Image performCut2(Image& dest, Image& psrc, int dx, int dy, float& cost){

	 vector<cut_node> cut;
     Image mask(psrc.width(),psrc.height());

     int drange = 10;

     // Determine overlapping area
	for (int x=0; x<psrc.width(); x++) for (int y=0; y<psrc.height(); y++)
        if (x+dx<dest.width() && y+dy<dest.height() && dest(x+dx,y+dy)!=0.)
            cut.push_back( cut_node (node_t(x,y),true) );

     if (cut.size()==0){
         cost = 0.;	return mask;
     }

    for (unsigned int k=0; k<cut.size(); k++){
        node_t p = cut[k].first;

		bool fg=false;
        if ((p.y==0 || p.x==0 || p.y==psrc.height()-1 || p.x==psrc.width()-1)){
			int nx; int ny;

            nx = p.x; ny = p.y-1;
			if ( ny<0 && ny+dy>=0 && dest(dx+nx, dy+ny)!=0. )	fg = true;

			nx = p.x; ny = p.y+1;
			if ( ny>=psrc.height() && ny+dy<dest.height() && dest(dx+nx, dy+ny)!=0. )	fg = true;

			nx = p.x-1; ny = p.y;
			if ( nx<0 && nx+dx>=0 && dest(dx+nx, dy+ny)!=0. )	fg = true;

			nx = p.x+1; ny = p.y;
			if ( nx>=psrc.width() && nx+dx<dest.width() && dest(dx+nx, dy+ny)!=0. )	fg = true;

        }
        cost=0.;
        if (fg){ //Coming from sink
            mask(p.x,p.y) = vsSINK;
            cost+= (psrc(p.x,p.y)-dest.atXY(dx+p.x,dy+p.y))*(psrc(p.x,p.y)-dest.atXY(dx+p.x,dy+p.y));
            //for (int m=0; m<drange; m++)  if (mask.inBounds(p.x+m,p.y+m)) mask(p.x+m,p.y+m) = vSINK; //to comment
            //if (mask.inBounds(p.x+drange,p.y+drange))         mask(p.x+drange,p.y+drange) = vsSINK; //to comment
        }


	 }

    return mask;
}

Image performCut3(Image& dest, Image& psrc, int dx, int dy, float& cost){

	 //vector<cut_node> cut;
     Image mask(psrc.width(),psrc.height());
     cost=0.;

    Image tmp(mask);
    //if (cut.size()>0){

		for (int i=0; i<mask.width(); i++) for (int j=0; j<mask.height(); j++){
		     if (i+dx<dest.width() && j+dy<dest.height() && dest(i+dx,j+dy)!=0.){
                tmp(i,j) = vSINK;
			 }
			 else 	tmp(i,j) = vSOURCE;
        }

        for (int i=0; i<mask.width(); i++) for (int j=0; j<mask.height(); j++)
            if (tmp(i,j)==vSINK){
                mask(i,j)=vSINK;
                 if ( tmp.atXY(i-1,j)==vSOURCE || tmp.atXY(i+1,j)==vSOURCE || tmp.atXY(i,j-1)==vSOURCE || tmp.atXY(i,j+1)==vSOURCE ||
                      tmp.atXY(i-1,j-1)==vSOURCE || tmp.atXY(i+1,j-1)==vSOURCE || tmp.atXY(i-1,j+1)==vSOURCE ||tmp.atXY(i+1,j+1)==vSOURCE )
                        mask(i,j) = vsSINK;
            }
            else if (tmp(i,j)==vSOURCE){
                mask(i,j)=vSOURCE;
                if ( tmp.atXY(i-1,j)==vSINK || tmp.atXY(i+1,j)==vSINK || tmp.atXY(i,j-1)==vSINK || tmp.atXY(i,j+1)==vSINK ||
                      tmp.atXY(i-1,j-1)==vSINK || tmp.atXY(i+1,j-1)==vSINK || tmp.atXY(i-1,j+1)==vSINK ||tmp.atXY(i+1,j+1)==vSINK )
                        mask(i,j) = vsSOURCE;
            }

    //}

    return mask;
}
