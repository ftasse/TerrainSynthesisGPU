#include "dtts_patchsynthesis.h"
#include <stdlib.h>
#include <iostream>
#include <ctime>
#include <iomanip>


#include "dtts_match.h"

void normalizeTerrain( float * h_data, int w, int h )
{
	float min, max;
	float height;
	int i;

	min = h_data[0];
	max = h_data[0];

	//find the min/max values of the height temp_buffer
	for( i=1; i<w*h; i++ )
	{
		if( h_data[i]>max )
			max= h_data[i];

		else if( h_data[i]<min )
			min= h_data[i];
	}

	//find the range of the altitude
	if( max <= min )
		return;

	height= max-min;

	//scale the values to a range of 0-255 (because I like things that way)
	for( i=0; i<w*h; i++ )
		h_data[i]= ( ( h_data[i]-min )/height )*255.0f;
}

int main(int argc, char** argv)
{
	    initrand();
            clock_t tstart=0, tstop=0;
            tstart = clock();

            char* exemplar = NULL;
            char* sketchf = NULL;
            char* output = NULL;
            char ridges = 'r';
            int bsize=100;

             initrand();


             /*Terrain dem;
             Tree dem_features;
             dem.loadTerragen(argv[1]);
             int who = atoi(argv[2]);
             if (who)
                dem_features.runPPA(dem,'r',9);
             else
                dem_features.runPPA_bis(dem,'r',9);
             return 0;*/


            /*{
                int mindelta = 0;
                int maxdelta = 255;
                int w=200;
                int l=200;
                int size = l;
                Terrain A(w,l);
                int iterations=1000;
                cout<<"hi\n";

            float * temp_buffer;
            int i,x,z;
            float height, height_reducer;
            int rect_size = size;
            int ni, nj;
            int mi, mj;
            int pmi, pmj;
            float temp;
            //height_map_size = size;


            int roughness = 1;


            //allocate the memory for our height data
            //height_data.data= new unsigned char [height_map_size*height_map_size];
            temp_buffer= A.getPixels();

            if( roughness<0 )
                roughness*= -1;

            height		  = ( float )rect_size/2;
            height_reducer= ( float )powf(2, -1*roughness);



            //set the first value in the height field
            temp_buffer[0]= 0.0f;

            int cur = 1;

            //being the displacement process
            while( rect_size>0 )
            {
                for( int i=0; i<size; i+=rect_size )
                {
                    for( int j=0; j<size; j+=rect_size )
                    {
                        ni= ( i+rect_size )%size;
                        nj= ( j+rect_size )%size;

                        mi= ( i+rect_size/2 );
                        mj= ( j+rect_size/2 );

                        temp = randfloat( -height/2, height/2 );
                        temp_buffer[mi+mj*size]= ( float )( ( temp_buffer[i+j*size] + temp_buffer[ni+j*size] + temp_buffer[i+nj*size] + temp_buffer[ni+nj*size] )/4 + randfloat( -height/2, height/2 ) );
                    }
                }

                for( int i=0; i<size; i+=rect_size )
                {
                    for( int j=0; j<size; j+=rect_size )
                    {

                        ni= (i+rect_size)%size;
                        nj= (j+rect_size)%size;

                        mi= (i+rect_size/2);
                        mj= (j+rect_size/2);

                        pmi= (i-rect_size/2+size)%size;
                        pmj= (j-rect_size/2+size)%size;

                        //Calculate the square value for the top side of the rectangle
                        temp_buffer[mi+j*size]= ( float )( ( temp_buffer[i+j*size]	  +
                                                                temp_buffer[ni+j*size]	  +
                                                                temp_buffer[mi+pmj*size]	  +
                                                                temp_buffer[mi+mj*size] )/4+
                                                                randfloat( -height/2, height/2 ) );

                        //Calculate the square value for the left side of the rectangle
                        temp_buffer[i+mj*size]= ( float )( ( temp_buffer[i+j*size]	  +
                                                                temp_buffer[i+nj*size]	  +
                                                                temp_buffer[pmi+mj*size]	  +
                                                                temp_buffer[mi+mj*size] )/4+
                                                                randfloat( -height/2, height/2 ) );
                    }
                }

                //reduce the rectangle size by two to prepare for the next
                //displacement stage
                rect_size/= 2;

                //reduce the height by the height reducer
                height*= height_reducer;

                //if ( cur==1 || cur==2 || cur==3 || cur==4 || cur==10 || cur==100 || cur==500 || cur==1000)
                {
                    Terrain B(w,l);
                    for( z=0; z<size; z++ )
                {
                    for( x=0; x<size; x++ )
                    {
                        B(x,z)=A(x,z);
                    }
                }

                 normalizeTerrain(B.getPixels(),w,l);
                 cout<<cur<<endl;
                 B.savePGM("/tmp/fault.pgm");
                 cin.get();
                }


                cur++;
            }*/

            /*
            Fault formation
                float* temp_buffer;
                int current_iteration;
                int height;
                int rand_x1, rand_z1;
                int rand_x2, rand_z2;
                int dir_x1, dir_z1;
                int dir_x2, dir_z2;
                int x, z;
                int i;

                //allocate the memory for our height data
                //height_data.data = new unsigned char [size*size];
                temp_buffer= A.getPixels();

            for( current_iteration=0; current_iteration<iterations; current_iteration++ )
            {
                //calculate the height range (linear interpolation from maxdelta to
                //mindelta) for this fault-pass
                height= maxdelta - ( ( maxdelta-mindelta )*current_iteration )/iterations;

                //pick two points at random from the entire height map
                rand_x1= rand( )%size;
                rand_z1= rand( )%size;

                //check to make sure that the points are not the same
                do
                {
                    rand_x2= rand( )%size;
                    rand_z2= rand( )%size;
                } while ( rand_x2==rand_x1 && rand_z2==rand_z1 );


                //dir_x1, dir_z1 is a vector going the same direction as the line
                dir_x1= rand_x2-rand_x1;
                dir_z1= rand_z2-rand_z1;

                for( z=0; z<size; z++ )
                {
                    for( x=0; x<size; x++ )
                    {
                        //dir_x2, dir_z2 is a vector from rand_x1, rand_z1 to the current point (in the loop)
                        dir_x2= x-rand_x1;
                        dir_z2= z-rand_z1;

                        //if the result of ( dir_x2*dir_z1 - dir_x1*dir_z2 ) is "up" (above 0),
                        //then raise this point by height
                        if( ( dir_x2*dir_z1 - dir_x1*dir_z2 )>0 )
                            temp_buffer[( z*size )+x]+= ( float )height;
                    }
                }


                //erode terrain
                int cur = current_iteration+1;
                if ( cur==1 || cur==2 || cur==3 || cur==4 || cur==10 || cur==100 || cur==500 || cur==1000)
                {
                    Terrain B(w,l);
                    for( z=0; z<size; z++ )
                {
                    for( x=0; x<size; x++ )
                    {
                        B(x,z)=A(x,z);
                    }
                }

                 normalizeTerrain(B.getPixels(),w,l);
                 cout<<current_iteration<<endl;
                 B.savePGM("/tmp/fault.pgm");
                 cin.get();
                }


                //filterHeightField( temp_buffer, filter );
            }
        A.savePGM("/tmp/fault.pgm");

        return 0;
    }*/

    /*{
        Terrain s;
        Terrain d(150,150);
        s.loadPGM("/home/flora/Dropbox/terrain_data/sdem2.pgm");

        ofstream out("/home/flora/Desktop/data.in");
        out<<"P5 10 10 255"<<endl;

        for (int j=0; j<s.width(); j+=s.width()/16)
            for (int i=0; i<s.height(); i+=s.width()/16){

                out<<s(i,j)<<endl;
            }
        out.close();
        return 0;

        int size = s.width()/4;
        int osize = 3*(size/4);
        //cout<<size<<osize<<endl; cin.get();

        patch_synthesis(d, s,size,osize);

        tstop = clock();
        // time in (ms)
        std::cerr<< "Total elapsed CPU time: "<<(float)(tstop-tstart)/(float)(CLOCKS_PER_SEC)<< " (s)" << std::endl;

        d.savePGM("/home/flora/Dropbox/terrain_data/sdem_patch_based.pgm");
        return 0;
    }*/

   /*{
        Terrain s;
        //s.loadPGM("/home/flora/Dropbox/terrain_data/dem.pgm");
        s.loadTerragen("/home/flora/Dropbox/terrain_data/CD.ter");
        Image cand1 = s.get_crop(10,10,109,109);    cand1.savePGM("/tmp/cand1.pgm");
        Image cand2 = s.get_crop(200,160,299,259);   cand2.savePGM("/tmp/cand2.pgm");

        for (int i=0; i<100; i++) for (int j=0; j<100; j++) cand1(i,j)=255;
        for (int i=0; i<100; i++) for (int j=0; j<100; j++) cand2(i,j)=0.01;

        Terrain res(150,150);
        res.mscale = s.mscale;

        for (int i=0; i<100; i++) for (int j=0; j<100; j++)  res(i,j) = cand1(i,j);
        res.savePGM("/tmp/candres.pgm");

        Image tmp (150,150);
        for (int i=0; i<100; i++) for (int j=0; j<100; j++)  tmp(i,j) = cand1(i,j);
        for (int i=50; i<100; i++) for (int j=0; j<100; j++)  tmp(i,j) = (cand2(i-50,j));
        for (int i=100; i<150; i++) for (int j=0; j<100; j++)  tmp(i,j) = cand2(i-50,j);
        tmp.savePGM("/tmp/candtmp.pgm");

        patch_merging(&res,&cand2,50,0);

        res.maxval=s.maxval;

        int ti = 50;
        for (int j=1; j<100; j++){
            //cout<<res(ti,j)-255<<" "<<res(ti,j)<<endl;
        }

        //shepard(&res,&cand2,50,0,1,100/4);
        //poisson(&res,&cand2,50,0);
        //paste_cut(&res,&cand2,50,0);
        res.savePGM("/tmp/candresf.pgm");
        res.saveTerragen("/tmp/candresf.ter");

        return 0;

    }*/

    /*cost_t *tmp = new cost_t [4096];
    for (int k=0; k<4096; k++){
    	tmp[k].org.x = rand()%100;
    	tmp[k].org.y = rand()%100;
    	tmp[k].cost = (rand()%100)/100.0f;
    	cout<<tmp[k].org.x<<" "<<tmp[k].org.y<<" "<<tmp[k].cost<<endl;
    }
    cin.get();
    // = {cost(node_t(0,0),103.0), cost(node_t(0,0),103.0), cost(node_t(0,0),13.4), cost(node_t(0,0),10003.1), cost(node_t(0,0),95.6)};
    sort_host(tmp,4096);
    cout<<"Done!\n"; cin.get();
    for (int k=0; k<4096; k++)
    	cout<<tmp[k].org.x<<" "<<tmp[k].org.y<<" "<<tmp[k].cost<<endl;
    delete [] tmp;
    return 0;*/


    if (argc>3 && (argv[1][1]=='s' || argv[1][1]=='S') )
    {
        exemplar = argv[2];
        output  =  argv[3];

        int nw = -1, nh = -1;

        if (argc>5)
        {
            nw = atoi(argv[4]);
            nh = atoi(argv[5]);
        }

        if (argc>6)
            bsize = atoi(argv[6]);
        int osize = bsize/4;
        if (argc>7)
            osize = atoi(argv[7]);

        terrain_synthesis(exemplar,output,nw,nh,bsize,osize);
    }
    else if (argc<6)
    {

        exemplar = (char*) "/home/flora/Dropbox/terrain_data/grand_canyon-e.ter"; //grand_canyon-e
        sketchf  = (char*) "/home/flora/Dropbox/terrain_data/exm.pgm";
        output = (char*) "/home/flora/Dropbox/terrain_data/terrsyn.ter";
        bsize=80;
        ridges = 'v';
        terrain_synthesis(exemplar,sketchf,output,ridges,bsize);
    }
    else
    {
        ridges = argv[1][1];
        exemplar = argv[2];
        sketchf  =  argv[3];
        output = argv[4];
        bsize=atoi(argv[5]);
        terrain_synthesis(exemplar,sketchf,output,ridges,bsize);
    }

        // the stop time counter
    tstop = clock();
    // time in (ms)
    std::cerr<< "Total elapsed CPU time: "<<(float)(tstop-tstart)/(float)(CLOCKS_PER_SEC)<< " (s)" << std::endl;

    return 0;

    }
