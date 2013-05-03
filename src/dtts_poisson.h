#include "external/linear_solver/ArrayVector.h"
#include "external/linear_solver/SparseMatrix.h"
#include "external/linear_solver/LinearSolver.h"

//#include "dtts_image.h"
//using namespace Dtts;


typedef pair<Image,Image> Gradient;

Gradient get_gradient(Image& src){

    /*Image res1 = src.convolute(sobel_filterX,3);
    Image res2 = src.convolute(sobel_filterY,3);
    Gradient grad(res1,res2);
    return grad;*/
    Gradient grad(Image(src.width(),src.height()),Image(src.width(),src.height()));

    for (int i=0; i<src.width(); i++) for (int j=0; j<src.height(); j++){
        grad.first(i,j)  = src.atXY(i+1,j)-src.atXY(i,j);
        grad.second(i,j) = src.atXY(i,j+1)-src.atXY(i,j);
    }

    return grad;//Gradient(res1,res2);
}

Gradient get_gradient_r(Image& src){

    //Image res1 = src.convolute(sobel_filterX,3);
    //Image res2 = src.convolute(sobel_filterY,3);
    Gradient grad(Image(src.width(),src.height()),Image(src.width(),src.height()));

    for (int i=0; i<src.width(); i++) for (int j=0; j<src.height(); j++){
        grad.first(i,j)  = src.atXY(i,j)-src.atXY(i-1,j);
        grad.second(i,j) = src.atXY(i,j)-src.atXY(i,j-1);
    }

    return grad;//Gradient(res1,res2);
}

Image get_divergent(Gradient grad){
    Image div(grad.first.width(),grad.first.height());
    Image gradXX = get_gradient_r(grad.first).first;
    Image gradYY = get_gradient_r(grad.second).second;

    for (int i=0; i<grad.first.width(); i++) for (int j=0; j<grad.first.height(); j++){
        div(i,j) = gradXX(i,j)+gradYY(i,j);
    }
    return div;
}

Image get_divergent(Image& src){
    return get_divergent(get_gradient(src));
}

/*Image wire_laplacian(Image& nsrc, Image& ndest, Image& mask){
    Image nsrc_grad = get_divergent(nsrc);
    Image ndest_div =
}*/

void poisson_solver_fast(Image& dest, Image& psrc, Image& mask, int drange, int dx, int dy){
    int nwidth = psrc.width() + 2*drange;
    int nheight = psrc.height() + 2*drange;

    for (int x=0; x<psrc.width(); x++) for (int y=0; y<psrc.height(); y++){
        if (x+dx<dest.width() && y+dy<dest.height()){
             if (mask(x,y)<=vsSOURCE)    dest(x+dx,y+dy) = psrc(x,y);
        }
    }

    Image nsrc(nwidth,nheight);
    Image nmask(nwidth,nheight);

    for (int x=0; x<nsrc.width(); x++) for (int y=0; y<nsrc.height(); y++) {

            if ( x+dx-drange>=0 &&  x+dx-drange<dest.width() && y+dy-drange>=0 &&  y+dy-drange<dest.height())
                nsrc(x,y) = dest(x+dx-drange,y+dy-drange);

            if ( x-drange>=0 &&  x-drange<mask.width() && y-drange>=0 &&  y-drange<mask.height() && mask(x-drange,y-drange)!=0){
                nmask(x,y) = mask(x-drange,y-drange);
                //for (int j=-drange/2; j<=drange/2; j++)   for (int i=-drange/2; i<=drange/2; i++)
                //    if (nmask(x+i,y+j)==0.) nmask(x+i,y+j) = vPAD;
            }
        }

    int *pos = new int [nwidth*nheight];
    for (int k=0; k<nwidth*nheight; k++)  pos[k] = -1;
    uint N=0;

	for (int x=0; x<nmask.width(); x++) for (int y=0; y<nmask.height(); y++)
        if (nmask(x,y)!= 0. )
        {
            pos[x+y*nwidth] = N;
			N++;
		}

    if (N == 0) {
		cout << "Solver::solve: No masked pixels found (mask color is non-grey)\n";
		delete [] pos;
		return;
	}

    //

	  for (int x=0; x<nsrc.width(); x++) for (int y=0; y<nsrc.height(); y++) {
        if (nsrc(x,y)==0.)
                nsrc(x,y) = psrc.atXY(x-drange,y-drange);
        //if (nmask(x,y)==0.)       nmask(x,y) = vSINK;
     }

    //nmask.savePGM("data/external/nmask.pgm");

    Gradient grad = get_gradient(nsrc);

    for (int i=0; i<nmask.width(); i++) for (int j=0; j<nmask.height(); j++){
        if (nmask(i,j)==vsSINK) // nmask(i,j)==vsSOURCE)
        {   grad.first(i,j) = 0;//psrc.atXY(i-drange,j-drange) - psrc.atXY((i-1)-drange,j-drange); //0.;
            grad.second(i,j) = 0;//psrc.atXY(i-drange,j-drange) - psrc.atXY(i-drange,(j-1)-drange); //0.;
            //nmask(i,j)=0.;
	    }
    }

    /*nmask.savePGM("data/tmp/nmask.pgm");
    grad.first.savePGM("data/tmp/grad0.pgm");
    grad.second.savePGM("data/tmp/grad1.pgm");*/



   //Image div = grad[0].get_gradient(0,4)[0]+grad[1].get_gradient(0,4)[1];

    LinearSolver S;
	S.Init(N);
	int count = 0;

	//printf("Solve matrix %d\n", N);
    ///nsrc.savePGM("data/tmp/nsrc_bef.pgm");
    //grad0.savePGM("data/tmp/gradX.pgm");

    Image div = get_divergent(grad);
    //cout<<grad0.getMax()<<" "<<grad1.getMax()<<endl;

    for (int i=0; i<nmask.width(); i++) for (int j=0; j<nmask.height(); j++){
        if (nmask(i,j)==vsSINK  || nmask(i,j)==vsSOURCE)
        {   //div(i,j)=0.;
	    }
    }

	for (int x=0; x<nsrc.width(); x++) for (int y=0; y<nsrc.height(); y++)
        if (pos[x+y*nwidth]>-1) {

            S.b[count]=0.;

            if (y-1>=0 && pos[x+(y-1)*nwidth]>-1) {
				int colIndex = pos[x+(y-1)*nwidth];
		        S.PushElement(colIndex,count,-1.);
            }
			else{ // at the top boundary
					S.b[count] += (float) nsrc.atXY(x,y-1);
			}

			if (y+1<nheight && pos[x+(y+1)*nwidth]>-1) {
				int colIndex = pos[x+(y+1)*nwidth];
		        S.PushElement(colIndex,count,-1.);
            }
			else{ // at the top boundary
					S.b[count] += (float) nsrc.atXY(x,y+1);
			}

			if (x-1>=0 && pos[(x-1)+y*nwidth]>-1) {
				int colIndex = pos[(x-1)+y*nwidth];
		        S.PushElement(colIndex,count,-1.);
            }
			else{ // at the top boundary
					S.b[count] += (float) nsrc.atXY(x-1,y);
			}

			if (x+1<nwidth-1 && pos[(x+1)+y*nwidth]>-1) {
				int colIndex = pos[(x+1)+y*nwidth];
		        S.PushElement(colIndex,count,-1.);
            }
			else{ // at the top boundary
					S.b[count] += (float) nsrc.atXY(x+1,y);
			}

			S.PushElement(count, count,4.);
			//if (nmask(x,y)!=vsSINK)
            //S.b[count] -= (float) 4*nsrc(x,y) - (nsrc.atXY(x-1,y)+nsrc.atXY(x+1,y)+nsrc.atXY(x,y+1)+nsrc.atXY(x,y-1));
            //if (!(nmask(x,y)>0 && nmask(x,y)<vPAD))
            S.b[count] -= div(x,y); //(grad0.atXY(x,y)-grad0.atXY(x+1,y))+(grad1.atXY(x,y)-grad1.atXY(x,y+1));
			//cout<<x<<" "<<y<<" "<<S.b[count]<<endl;
            count++;

		}


	S.BiCGradSolve(300);

	for (int x=0; x<nsrc.width(); x++) for (int y=0; y<nsrc.height(); y++)
        if (pos[x+y*nwidth]>-1){
            int index=pos[x+y*nwidth];
            float nval =  (float) S.x[index];
		    nsrc(x,y) =nval;
		    nmask(x,y)=vPAD;
        }
    ///nsrc.savePGM("data/tmp/nsrc_aft.pgm");
    ///nmask.savePGM("data/tmp/nmask.pgm");

    cout<<"Size of unknowns: "<<N<<endl;

    for (int x=0; x<nsrc.width(); x++) for (int y=0; y<nsrc.height(); y++)
            if ( x+dx-drange>=0 &&  x+dx-drange<dest.width() && y+dy-drange>=0 &&  y+dy-drange<dest.height() && dest(x+dx-drange,y+dy-drange)!=0.)
                dest(x+dx-drange,y+dy-drange)=nsrc(x,y);
    delete [] pos;
}

void poisson_solver_slow(Image& dest, Image& psrc, Image& mask, int drange, int dx, int dy){
    int nwidth = psrc.width() + 2*drange;
    int nheight = psrc.height() + 2*drange;

    /*for (int x=0; x<psrc.width(); x++) for (int y=0; y<psrc.height(); y++){
        if (x+dx<dest.width() && y+dy<dest.height()){
             if (mask(x,y)<=vsSOURCE)    dest(x+dx,y+dy) = psrc(x,y);
        }
    }*/

    Image nsrc(nwidth,nheight);
    Image nmask(nwidth,nheight);
    Image ntmp(nwidth,nheight);

    for (int x=0; x<nsrc.width(); x++) for (int y=0; y<nsrc.height(); y++) {

            if ( x+dx-drange>=0 &&  x+dx-drange<dest.width() && y+dy-drange>=0 &&  y+dy-drange<dest.height())
                nsrc(x,y) = dest(x+dx-drange,y+dy-drange);

            ntmp(x,y) = psrc.atXY(x-drange,y-drange);

            if ( x-drange>=0 &&  x-drange<mask.width() && y-drange>=0 &&  y-drange<mask.height()){// && mask(x-drange,y-drange)!=0){
                nmask(x,y) = mask(x-drange,y-drange);
                if (nmask(x,y)<=vsSINK)
                    nmask(x,y)=vPAD;    //nmask(x,y)>0 &&
                //for (int j=-drange/2; j<=drange/2; j++)   for (int i=-drange/2; i<=drange/2; i++)
                //    if (nmask(x+i,y+j)==0.) nmask(x+i,y+j) = vPAD;
            }
        }

    int *pos = new int [nwidth*nheight];
    for (int k=0; k<nwidth*nheight; k++)  pos[k] = -1;
    uint N=0;

	for (int x=0; x<nmask.width(); x++) for (int y=0; y<nmask.height(); y++)
        //if (nmask(x,y)!= 0. )
        if (nmask(x,y) == vPAD )
        {
            pos[x+y*nwidth] = N;
			N++;
		}

    if (N == 0) {
		//cout << "Solver::solve: No masked pixels found (mask color is non-grey)\n";
		delete [] pos;
		return;
	}

    //nmask.savePGM("data/tmp/nmask.pgm");

	  for (int x=0; x<nsrc.width(); x++) for (int y=0; y<nsrc.height(); y++) {
        if (nsrc(x,y)==0.)
                nsrc(x,y) = psrc.atXY(x-drange,y-drange);
        //if (nmask(x,y)==0.)       nmask(x,y) = vSINK;
     }

    ///ntmp.savePGM("data/tmp/ntmp.pgm");

    Image div (nwidth,nheight);//= get_divergent(get_gradient(ntmp));

    for (int i=0; i<nmask.width(); i++) for (int j=0; j<nmask.height(); j++){
        //div(i,j)= 4*nsrc.atXY(i,j)-(nsrc.atXY(i,j-1)+nsrc.atXY(i,j+1)+nsrc.atXY(i-1,j)+nsrc.atXY(i+1,j));
        //if (nmask(i,j)>0 && nmask(i,j)<vPAD)
        //if (nmask(i,j)==vPAD)// || nmask(i,j)==vsSOURCE)
        {
            div(i,j)= 4*ntmp.atXY(i,j)-(ntmp.atXY(i,j-1)+ntmp.atXY(i,j+1)+ntmp.atXY(i-1,j)+ntmp.atXY(i+1,j));
            //nmask(i,j)=0.;
	    }
    }

   ///nmask.savePGM("data/tmp/nmask.pgm");
   //Image div = grad[0].get_gradient(0,4)[0]+grad[1].get_gradient(0,4)[1];

   LinearSolver S;
	S.Init(N);
	int count = 0;

	//printf("Solve matrix %d\n", N);
    ///nsrc.savePGM("data/tmp/nsrc_bef.pgm");
    //grad0.savePGM("data/tmp/gradX.pgm");

	for (int x=0; x<nsrc.width(); x++) for (int y=0; y<nsrc.height(); y++)
        if (pos[x+y*nwidth]>-1) {

            S.b[count]=0.;

            if (y-1>=0 && pos[x+(y-1)*nwidth]>-1) {
				int colIndex = pos[x+(y-1)*nwidth];
		        S.PushElement(colIndex,count,-1.);
            }
			else{ // at the top boundary
					S.b[count] += (float) nsrc.atXY(x,y-1);
			}

			if (y+1<nheight && pos[x+(y+1)*nwidth]>-1) {
				int colIndex = pos[x+(y+1)*nwidth];
		        S.PushElement(colIndex,count,-1.);
            }
			else{ // at the top boundary
					S.b[count] += (float) nsrc.atXY(x,y+1);
			}

			if (x-1>=0 && pos[(x-1)+y*nwidth]>-1) {
				int colIndex = pos[(x-1)+y*nwidth];
		        S.PushElement(colIndex,count,-1.);
            }
			else{ // at the top boundary
					S.b[count] += (float) nsrc.atXY(x-1,y);
			}

			if (x+1<nwidth-1 && pos[(x+1)+y*nwidth]>-1) {
				int colIndex = pos[(x+1)+y*nwidth];
		        S.PushElement(colIndex,count,-1.);
            }
			else{ // at the top boundary
					S.b[count] += (float) nsrc.atXY(x+1,y);
			}

			S.PushElement(count, count,4.);
			//if (nmask(x,y)!=vsSINK)
            //S.b[count] -= (float) 4*nsrc(x,y) - (nsrc.atXY(x-1,y)+nsrc.atXY(x+1,y)+nsrc.atXY(x,y+1)+nsrc.atXY(x,y-1));
            //if (!(nmask(x,y)>0 && nmask(x,y)<vPAD))
            S.b[count] += div(x,y);//(grad0.atXY(x,y)-grad0.atXY(x+1,y))+(grad1.atXY(x,y)-grad1.atXY(x,y+1));
			//cout<<x<<" "<<y<<" "<<S.b[count]<<endl;
            count++;

		}


	S.BiCGradSolve(300);

	for (int x=0; x<nsrc.width(); x++) for (int y=0; y<nsrc.height(); y++)
        if (pos[x+y*nwidth]>-1){
            int index=pos[x+y*nwidth];
            float nval =  (float) S.x[index];
		    nsrc(x,y) =nval;
        }
    ///nsrc.savePGM("data/tmp/nsrc_aft.pgm");

    //cout<<"Size of unknowns: "<<N<<endl;
    for (int x=0; x<psrc.width(); x++) for (int y=0; y<psrc.height(); y++){
        if (x+dx<dest.width() && y+dy<dest.height()){
             if (mask(x,y)<=vsSOURCE)    dest(x+dx,y+dy) = psrc(x,y);
        }
    }

    for (int x=0; x<nsrc.width(); x++) for (int y=0; y<nsrc.height(); y++)
            if ( x+dx-drange>=0 &&  x+dx-drange<dest.width() && y+dy-drange>=0 &&  y+dy-drange<dest.height() && dest(x+dx-drange,y+dy-drange)!=0.)
                dest(x+dx-drange,y+dy-drange)=nsrc(x,y);
    delete [] pos;
}

