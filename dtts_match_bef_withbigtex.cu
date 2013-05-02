
__host__ __device__ float get_diff(float* dest, float* bigtex, int bsize,  int cx, node_t p, node_t q)
{
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
        if (steep)
        {
	    int ide = y + x*bsize;
            float val = dest[ide]-bigtex[cx+(ide*bigw)];
            sum+= val*val;
        }
        else
        {
	    int ide = x + y*bsize;
            float val = dest[ide]-bigtex[cx+(ide*bigw)];
   
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

    return sqrt(sum)/count;
}
__host__ __device__ float get_diff_profile(float* dest, float* bigtex, int bsize, int cx, node_t* cand_br1, int csize)
{
    float total=0.;
    int mw=bsize, mh=bsize;
    node_t mid(mw/2,mh/2);

    for (int k=0; k<csize; k++)
    {
        total += get_diff(dest,bigtex,bsize,cx,cand_br1[k],mid);
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
        total += get_diff(dest,bigtex,bsize,cx,qnode1, mid);
        total += get_diff(dest,bigtex,bsize,cx,mid, qnode2);
    }

    return total;
}

__host__ __device__ float ssdf(float* dest, float* bigtex, int bsize, int cx)
{
    float sum = 0.;
    int count = 0; 
    //int ide = 0;
    //int alln = bigwidth/(bsize*bsize);
    for (int x=0; x<bsize*bsize; x++)
            if ( dest[x]!=0.)
            {
                //count++;
                float val = dest[x]-bigtex[cx+x*bigw];
                sum+= val*val;
                count++;
            }
    count++;
    return sqrt(sum)/count;
}

void getCand(Image& cand, float* bigtex, int cx)
{
		
    for (int j=0; j<cand.height(); j++)
        for (int i=0; i<cand.width(); i++){
            //cout<<cx<<" "<<bigw<<" "<<bigtex[cx]<<endl;
            cand(i,j)=bigtex[cx];
            cx+=bigw;
            
            }
    
    //cand.savePGM("/tmp/cand_img.pgm"); cin.get();
}


Image findCand(Image& dest, float* bigtex, int bigwidth, cost_t* candidates, int csize, cost_t* prev, int bsize, int dx, int dy)
{
    vector<cost_t> choices;
    //sort(candidates,candidates+csize);
    int s = csize;
    if (s>KNUM) s=KNUM;

    for (int k=0; k<s; k++)
    {
        Image cand(bsize,bsize);
        getCand(cand,bigtex,candidates[k].tnode.x);
        candidates[k].cost = graphCut_cost(&dest,&cand,dx,dy);
        choices.push_back(candidates[k]);
    }

    sort(choices.begin(),choices.end());
    cost_t choice = choices[0];
    
    Image res(bsize,bsize);
    getCand(res,bigtex,choice.tnode.x);
    
    *prev = choice;

    return res;
}

Image findCand_dev(Image& dest, float* bigtex, thrust::device_vector<cost_t> candidates_dev, int csize, cost_t* prev, int bsize, int dx, int dy)
{
    int s = csize;
    if (s>KNUM) s=KNUM;
    
    cost_t choice;    float mini = 10*INF;

    for (int k=0; k<s; k++)
    {
        Image cand(bsize,bsize);
	cost_t cur = candidates_dev[k];
        getCand(cand,bigtex,cur.tnode.x);
        float cost = graphCut_cost(&dest,&cand,dx,dy);
	if (cost<mini){
		mini = cost;
		choice = cur;
	}
    }
    
    Image res(bsize,bsize);
    getCand(res,bigtex,choice.tnode.x);
    
    *prev = choice;

    return res;
}

__host__ __device__ float getCost_noFeature(float* dest, float* bigtex, float* dem_vars, float* usr_var, cost_t dem, int bsize)
{
    int cx = dem.tnode.x;
    //int cy = dem.tnode.y;
    float tmp = 10*ssdf(dest,bigtex,bsize,cx);
    if (use_noisestat) tmp+=0.0001*compare_variances(usr_var,dem_vars,dem.vpos,NLEVEL);
    return tmp;
}

__global__ void ComputeCosts_noFeature_Kernel( float* dest, float* bigtex, cost_t* candidates, int csize, float* dem_vars, float* usr_var, cost_t* prev, int bsize){
	//__shared__ float* dest = destg;
	int k = blockDim.x * blockIdx.x + threadIdx.x;
	if (k<csize){
		if ( (!candidates[k].skip) && k>5)
		{
		    candidates[k].cost = getCost_noFeature(dest, bigtex, dem_vars, usr_var, candidates[k], bsize);
		}
		else
		{
		    candidates[k].cost = INF;
		}
	}
	
}

__global__ void buildCandidates_noFeature_Kernel(cost_t* candidates, float* bigtex, float* src_ptr, node_t* dnodes,int nsize, int src_w, int src_h, int bsize){
    
    int k = blockDim.x * blockIdx.x + threadIdx.x;
    
   if (k<nsize*(360/DROT+ DMIR))
    {
    
    int kn = k%(nsize); 
    int rt = k/(nsize); 
    
    bigw = nsize*(360/DROT+ DMIR);
    
    
        node_t pnode = dnodes[kn];
        int rx = pnode.x;
        int ry = pnode.y;
        
        int mid = bsize/2;
        int rs = (360/DROT+ DMIR);
        int rot = rt*DROT;

       	    cost_t c(node_t(rx,ry),0,rot,0,0);
        
            {
            	int kpos = (kn*(rs))+rt;
                c.vpos = kpos*3;

                int cx = kpos;  //for coleasced memory access
                int cy = 0;
                c.tnode = node_t(cx,cy);

                c.skip = false;
		for (int j=0; j<bsize; j++)
                    for (int i=0; i<bsize; i++)
                    {
                        //cout<<cx+i<<"/"<<cy+i<<": "<<cand(i,j)<<endl;
                        int ri = 0;
                        int rj = 0;
                        float candv=0;
			if (rt<360/DROT){
				 float ni=0;
                        	 float nj=0;
		                ni = (rx+mid) + ((i - mid)*cos_int(rot)) + ((j - mid) * sin_int(rot));
		                nj = (ry+mid) - ((i - mid)*sin_int(rot)) + ((j - mid) * cos_int(rot));
		                if (ni>=0 && nj>=0 && ni<src_w && nj<src_h)
		              		candv = cubicInterpol(src_ptr,src_w,src_h,ni,nj);    
                            }
                        else if (rot==360/DROT){
                        	ri = (rx+bsize-i-1);
                            	rj = (ry+j);
                            	 if( ri>=0 && rj>=0 && ri<src_w && rj<src_h)
                         		candv = src_ptr[ ri + rj*src_w  ];
                        }
                        else {
                            ri = (rx+i);
                            rj = (ry+bsize-j-1);
                             if( ri>=0 && rj>=0 && ri<src_w && rj<src_h)
                         		candv = src_ptr[ ri + rj*src_w  ];
                        
                        }                            
                        
                        (bigtex)[cx] = candv;
                        cx+= (rs*nsize); //for coleasced memory access
                        if (candv<0.0001){
                        	c.skip = false;
                        	candidates[kpos] = (c);
                        	//return;                        	
                        }
                    }

                candidates[kpos] = (c);
        }
    }
}

void match_noFeature_bef(Terrain& dest, Image& src, Image& target, node_list dem_nodes, vector<Image>& tar_pyr, vector<Image>&  src_pyr, int bsize, int osize)
{
    clock_t start_t, end_t, s_tmp, e_tmp;
    start_t = clock();
    
    match_time = 0;
    paste_time = 0;
    get_target = 0;


    int nsize = dem_nodes.size();

   
    
    {
    end_t = clock();
    float elapsed = ((float)( end_t - start_t )) /CLOCKS_PER_SEC;
    cout<<"Get pyramid elapsed time: "<<elapsed<<" s.\n";
    start_t = clock();
    }

   float* variances = get_noise_stats_Feature(src_pyr, dem_nodes, bsize);
   
   {
    end_t = clock();
    float elapsed = ((float)( end_t - start_t )) /CLOCKS_PER_SEC;
    cout<<"Get noise stats elapsed time: "<<elapsed<<" s.\n";
    start_t = clock();
    }

   node_t* dnodes = new node_t [nsize];
   int rs = (360/DROT);

    int bigwidth  = nsize*bsize;
    int bigheight = (rs+DMIR)*bsize;
    bigw = nsize*(rs+DMIR);
    Image big(bigwidth,bigheight);
    float* bigtex = big.getPixels();
    float* src_ptr = src.getPixels();

    cost_t* candidates = new cost_t[nsize*(rs+DMIR)];
    float* dem_vars = new float [nsize*(rs+DMIR)*3];

    {
	    int count=0;
	    int rv = 0;
	    for (node_list::const_iterator it = dem_nodes.begin(); it != dem_nodes.end(); it++ ){
	    	dnodes[count] = *it;
	    	
	    	for (int rot=0; rot<360; rot+=DROT)
            {
                for (int k=0; k<NLEVEL; k++)
                    dem_vars[rv++] = variances[count*NLEVEL+k];
            }
            for (int m=0; m<DMIR; m++){
                int mir = m+1;
                for (int k=0; k<NLEVEL; k++)
                    dem_vars[rv++] = variances[count*NLEVEL+k];
            }
            count++;
	    }
    }

    {
    end_t = clock();
    float elapsed = ((float)( end_t - start_t )) /CLOCKS_PER_SEC;
    cout<<"Elapsed time: "<<elapsed<<" s.\n";
    start_t = clock();
    }

    float cand_time = 0;
    cout<<"Prepare candidates! "<< nsize<<"\n";
    
    float* src_dev;		cudaMalloc((void**) &src_dev, sizeof(float)*src.getWidth()*src.getHeight());
    	
    float* bigtex_dev;		cudaMalloc((void**) &bigtex_dev, sizeof(float)*bigwidth*bigheight);
    float* dem_vars_dev;	cudaMalloc((void**) &dem_vars_dev, sizeof(float)*nsize*rs*3);
    
    //cost_t* candidates_dev;	cudaMalloc((void**) &candidates_dev, sizeof(cost_t)*nsize*rs);
    node_t* dnodes_dev;		cudaMalloc((void**) &dnodes_dev, sizeof(node_t)*nsize);
    
    cudaMemcpy(src_dev, src.getPixels(), sizeof(float)*src.getWidth()*src.getHeight(), cudaMemcpyHostToDevice);
    //cudaMemcpy(bigtex_dev, big.getPixels(), sizeof(float)*big.getWidth()*big.getHeight(), cudaMemcpyHostToDevice);
    
    
    cudaMemcpy(dem_vars_dev, dem_vars, sizeof(float)*nsize*(rs+DMIR)*3, cudaMemcpyHostToDevice);
    thrust::device_vector<cost_t> candidates_dev(nsize*(rs+DMIR));//cudaMemcpy(candidates_dev, candidates, sizeof(cost_t)*nsize*rs, cudaMemcpyHostToDevice);
    cudaMemcpy(dnodes_dev, dnodes, sizeof(node_t)*nsize, cudaMemcpyHostToDevice);
    
    
    //cout<<"Start matching non-feature\n"<<endl;
	
    {
    
    	dim3 dimBlock(32,rs+DMIR);
	dim3 dimGrid( (nsize / dimBlock.x)+1,1); //, A.height / dimBlock.y
	
	int threadsPerBlock = 256;
        int blocksPerGrid =  ((nsize*(rs+DMIR)) / threadsPerBlock)+1;
	
	 cout<<blocksPerGrid<<" "<<threadsPerBlock<<endl;
	
	//buildCandidates_noFeature_Kernel<<<dimGrid, dimBlock>>> (thrust::raw_pointer_cast(&candidates_dev[0]), bigtex_dev, src_dev, dnodes_dev,nsize, src.getWidth(),src.getHeight(),bsize);
	buildCandidates_noFeature_Kernel<<<blocksPerGrid, threadsPerBlock>>> (thrust::raw_pointer_cast(&candidates_dev[0]), bigtex_dev, src_dev, dnodes_dev,nsize, src.getWidth(),src.getHeight(),bsize);
	
		
    	//cudaMemcpy(candidates, candidates_dev, sizeof(cost_t)*nsize*rs, cudaMemcpyDeviceToHost);
    	//cudaMemcpy(big.getPixels(), bigtex_dev, sizeof(float)*bigwidth*bigheight, cudaMemcpyDeviceToHost);

    	end_t = clock();
    	float elapsed = ((float)( end_t - start_t )) /CLOCKS_PER_SEC;
    	cout<<"Elapsed time: "<<elapsed<<" s.\n";
    	cand_time += elapsed;
    	start_t = clock();
    }

    //big.savePGM("/tmp/bigtmp_noFeature.pgm");

    cost_t* prev = new cost_t(node_t(-1,-1),0,0,0,0);
    int dx = -5*bsize, dy=-5*bsize;

	int cnum = 0;
	for (int k=0; k<nsize*(rs+DMIR); k++)
		if (!candidates[k].skip)
			cnum++;
    cout<<"Start matching non-feature patches! from "<<cnum<<"\n";
    
    //return;
    
    float* usr_var_dev;    cudaMalloc((void**) &usr_var_dev,sizeof(float)*NLEVEL);
    cost_t* prev_dev; 	   cudaMalloc((void **) &prev_dev,sizeof(cost_t));
    float* ucand_dev; 	   cudaMalloc((void **) &ucand_dev, sizeof(float)*bsize*bsize);

    cnum = 0;
    
    dim3 dimBlock(BLOCK_SIZE,BLOCK_SIZE);
    dim3 dimGrid( (dest.width() / dimBlock.x)+1,(dest.height() / dimBlock.y)+1); //, A.height / dimBlock.y
    //cout<<dimGrid.x<<"/"<<dimGrid.y<<": "<<dimBlock.x<<"/"<<dimBlock.y<<endl;
    int ide = 220 + 74*dest.width();
    
    //cout<<ide<<" "<<cand_var[ide]<<endl;
    list<cost_t> omega;
    float* conf = new float [dest.width()*dest.height()];
    for (int x=0; x<dest.width(); x++)	for (int y=0; y<dest.height(); y++)	
    	if (dest(x,y)>BG)
    		conf[x+y*dest.width()] = 1;
    while (true)
    {
    	int  threadsPerBlock = 256;
    	int blocksPerGrid =  ((nsize*(rs+DMIR)) / threadsPerBlock)+1;
    	
        s_tmp = clock();
	//getNextTarget<<<dimGrid,dimBlock>>> (dest_dev,cand_var_dev,dest.width(), dest.height(), bsize, dx+bsize/2, dy+bsize/2);
	//cudaMemcpy(cand_var, cand_var_dev, sizeof(float)*dest.getWidth()*dest.getHeight(), cudaMemcpyDeviceToHost);
	
	getNextTarget2 (dest.getPixels(), conf, omega,dest.width(), dest.height(), bsize, dx+bsize/2, dy+bsize/2);
	if (omega.size()==0) break;
	
	//getNextTarget(dest.getPixels(), cand_var, dest.width(), dest.height(), bsize, dx+bsize/2, dy+bsize/2);
	//cout<<ide<<" "<<cand_var[ide]<<endl;	
	//cin.get();	
	
	int maxi = -1;
	int xmaxi = -1;
	int ymaxi = -1;
		float maxv = -INF, tmp; 
		
		for (list<cost_t>::iterator it=omega.begin(); it!=omega.end(); it++){
			tmp = (*it).cost;
			if (tmp>maxv){
				node_t p = (*it).org;
				maxv = tmp;
				xmaxi = p.x;
				ymaxi = p.y;
			}
		}
		

		
		//cout<<xmaxi<<"/"<<ymaxi<<" --> "<<maxv<<endl;
	e_tmp = clock();
	get_target+=mstimer(s_tmp,e_tmp);	

        //if (xmaxi<0) break;
        //dy=(maxi/dest.width())-bsize/2;
        //dx=(maxi%dest.width())-bsize/2;
        dx=(xmaxi)-bsize/2;
        dy=(ymaxi)-bsize/2;
        
        s_tmp = clock();

       vector<float> usr_var = noise_variances(tar_pyr,dx,dy,bsize);
       Image ucand = dest.get_crop(dx,dy,dx+bsize-1,dy+bsize-1);
       cudaMemcpy(usr_var_dev, &usr_var[0], sizeof(float)*NLEVEL, cudaMemcpyHostToDevice);
       
       cudaMemcpy(ucand_dev, ucand.getPixels(), bsize*bsize*sizeof(float), cudaMemcpyHostToDevice);

       ComputeCosts_noFeature_Kernel<<<blocksPerGrid, threadsPerBlock>>>(ucand_dev, bigtex_dev,thrust::raw_pointer_cast(&candidates_dev[0]) , nsize*(rs+DMIR), dem_vars_dev, usr_var_dev, prev_dev, bsize);
       thrust::sort(candidates_dev.begin(),candidates_dev.end(),comp);

                                        //cudaMemcpy(candidates, candidates_dev, sizeof(cost_t)*nsize*rs, cudaMemcpyDeviceToHost);
                                       //Image patch = findCand(dest,bigtex,bigwidth,candidates,nsize*rs, prev,bsize, dx,dy);
       Image patch = findCand_dev(dest,bigtex,candidates_dev,nsize*(rs+DMIR), prev,bsize,dx,dy);
       
        cudaMemcpy(prev_dev, prev, sizeof(cost_t), cudaMemcpyHostToDevice);
       
       e_tmp = clock();
       match_time+=mstimer(s_tmp,e_tmp);
       
       s_tmp = clock();
       paste_patch(dest,patch, bsize/10, dx, dy);
       //dest.savePGM("/tmp/res_tmp_gpu.pgm");
       //cudaMemcpy(dest_dev, dest.getPixels(), sizeof(float)*dest.width()*dest.height(), cudaMemcpyHostToDevice);
       
       e_tmp = clock();
       paste_time+=mstimer(s_tmp,e_tmp);
       
       		for (int i=0; i<bsize ;i++)	for (int j=0; j<bsize ;j++)	
			 if (xmaxi+i>=0 && ymaxi+j>=0 && xmaxi+i<dest.width() && ymaxi+j<dest.height() && dest(xmaxi+i,ymaxi+j)<=BG)
			 	conf[(xmaxi+i)+(ymaxi+j)*dest.width()] = conf[xmaxi+ymaxi*dest.width()];

        cnum++;
    }
    cout<<"Number of targets: "<<cnum<<endl;
    delete [] conf;

    {
    end_t = clock();
    float elapsed = ((float)( end_t - start_t )) /CLOCKS_PER_SEC;
    cout<<"Elapsed time: "<<elapsed<<" s.\n";
    start_t = clock();
    }
    
    int threadsPerBlock = 256;
    int blocksPerGrid =  ((nsize*(rs+DMIR)) / threadsPerBlock)+1;
    
    	bool havdata = false;
    	for (int id = 0; id<dest.height()*dest.width(); id++)
		if (dest.getPixels()[id]>BG){
			havdata = true;
			break;
		}
				

	if (cnum==0){
		
		cout<<"Number of targets: "<<cnum<<endl;
		while (true){
		bool finish = true;
		for (int y = 0; y<dest.height(); y++)
			for (int x = 0; x<dest.width(); x++)
				if (dest(x,y)<=BG && ( (!havdata) || onBoundary(dest.getPixels(),dest.width(),dest.height(),x,y)))
				{
					//dy=dx;
					
					finish = false;
					dx = x-bsize/2;
					dy = y-bsize/2;
        
				s_tmp = clock();

			       vector<float> usr_var = noise_variances(tar_pyr,dx,dy,bsize);
			       Image ucand = dest.get_crop(dx,dy,dx+bsize-1,dy+bsize-1);
			       cudaMemcpy(usr_var_dev, &usr_var[0], sizeof(float)*NLEVEL, cudaMemcpyHostToDevice);
			       cudaMemcpy(prev_dev, prev, sizeof(cost_t), cudaMemcpyHostToDevice);
			       cudaMemcpy(ucand_dev, ucand.getPixels(), bsize*bsize*sizeof(float), cudaMemcpyHostToDevice);

			       ComputeCosts_noFeature_Kernel<<<blocksPerGrid, threadsPerBlock>>>(ucand_dev, bigtex_dev,thrust::raw_pointer_cast(&candidates_dev[0]) , nsize*(rs+DMIR), dem_vars_dev, usr_var_dev, prev_dev, bsize);
			       thrust::sort(candidates_dev.begin(),candidates_dev.end(),comp);

						                //cudaMemcpy(candidates, candidates_dev, sizeof(cost_t)*nsize*rs, cudaMemcpyDeviceToHost);
						               //Image patch = findCand(dest,bigtex,bigwidth,candidates,nsize*rs, prev,bsize, dx,dy);
			       Image patch = findCand_dev(dest,bigtex,candidates_dev,nsize*(rs+DMIR), prev,bsize,dx,dy);
			        cudaMemcpy(prev_dev, prev, sizeof(cost_t), cudaMemcpyHostToDevice);
			       e_tmp = clock();
			       match_time+=mstimer(s_tmp,e_tmp);
			       
			       s_tmp = clock();
			       paste_patch(dest,patch, bsize/10, dx, dy);
			       //dest.savePGM("/tmp/res_tmp_gpu.pgm");
			       //cudaMemcpy(dest_dev, dest.getPixels(), sizeof(float)*dest.width()*dest.height(), cudaMemcpyHostToDevice);
			       
			       e_tmp = clock();
			       paste_time+=mstimer(s_tmp,e_tmp);

				cnum++;
					//dest.savePGM("/tmp/tmp.pgm");
				}
			
			    if (finish)	break;
			}

	

		end_t = clock();
		float elapsed = ((float)( end_t - start_t )) /CLOCKS_PER_SEC;
		cout<<"Elapsed time: "<<elapsed<<" s.\n";
	}


    cudaFree(usr_var_dev);
    cudaFree(prev_dev);
    delete prev;
    cudaFree(ucand_dev);

    end_t = clock();
    float elapsed = ((float)( end_t - start_t )) /CLOCKS_PER_SEC;
    cout<<"Elapsed time: "<<elapsed<<" s.\n";

    delete [] candidates;
    delete [] dnodes;
    delete [] dem_vars;
    delete [] variances;
    
    //cudaFree(candidates_dev);
    cudaFree(dnodes_dev);
    cudaFree(dem_vars_dev);
    cudaFree(src_dev);
    cudaFree(bigtex_dev);
    
    cerr<<"\n\n*********** Non Feature matching GPU*******************\n";
    print_times();
    cerr<<" Candidates set: "<<cand_time<<"s\n";
    cerr<<" Number  of targets: "<<cnum<<"\n";
    cerr<<" Number  of cands: "<<nsize*(rs+DMIR)<<"\n";
    cerr<<"*********** End *******************\n\n";
}

__host__ __device__ float getCost_Feature(float* dest, float* target, float* bigtex, float* dem_vars, node_t* dem_leafs, float* usr_var, node_t* uleafs, cost_t dem, int bsize)
{
    float tmp = 0.;

    int cx = dem.tnode.x;
    int cy = dem.tnode.y;


    //Image cand = bigtex.get_crop(dem.tnode.x,dem.tnode.y,dem.tnode.x+(bsize-1),dem.tnode.y+(bsize-1));
    //if (use_bend && usr.candleafs.size()>2) tmp+= 100*tps_img(cand,uleafs,dem_leafs,dem.lpos,bsize,dem.lsize);
    //if (use_cut) tmp+=graphCut_cost(&dest,&cand,dx,dy);
    //if (use_bend && dem.lsize>=3) tmp+= 1000*get_tps(uleafs,dem_leafs,dem.lpos,bsize,dem.lsize) ;
   
	if (use_noisestat)  tmp+=0.001*compare_variances(usr_var,dem_vars,dem.vpos,NLEVEL);
    
	if (use_angle) tmp+= 2*getDiffAng(uleafs,dem_leafs,dem.lpos,bsize,dem.lsize);
    //cout<<"hi\n";

	if (use_profile)  tmp+= 5*get_diff_profile(target,bigtex,bsize,cx,uleafs,dem.lsize);

    //cout<<"hey\n";
    
       if (use_ssd) tmp+= ssdf(dest,bigtex,bsize,cx);

    /*float t = ssdf(dest,bigtex,dsize,bsize,dx,dy,cx,cy);
    if (t>0){

        cout<<0.01*compare_variances(usr_var,dem_vars,dem.vpos,NLEVEL)<<endl;
        cout<<getDiffAng(uleafs,dem_leafs,dem.lpos,bsize,dem.lsize)<<endl;
        cout<<get_diff_profile(target,bigtex,bsize,dx,dy,cx,cy,uleafs,dem.lsize)<<endl;
        cout<<t/bsize<<endl;
    }*/

    return tmp;
}

__global__ void ComputeCosts_Feature_Kernel( float* dest, float* target, float* bigtex, cost_t* candidates, int csize, float* dem_vars, node_t* dem_leafs, float* usr_var, node_t* uleafs, cost_t* prev, int bsize, int lsize){

	int k = blockDim.x * blockIdx.x + threadIdx.x;
	if (k<csize){
		if ( (!candidates[k].skip) && lsize==candidates[k].lsize && k>5)
		{
		    candidates[k].cost = getCost_Feature(dest, target, bigtex, dem_vars, dem_leafs, usr_var, uleafs, candidates[k], bsize);
		}
		else
		{
		    candidates[k].cost = INF;
		}
	}
	
}

__global__ void buildCandidates_Feature_Kernel(cost_t* candidates, float* bigtex, float* src_ptr, node_t* dnodes, node_t* dem_lsizes, int nsize, int src_w, int src_h, int bsize){
    
   int k = blockDim.x * blockIdx.x + threadIdx.x;
    
   if (k<nsize*(360/DROT+ DMIR))
    {
    int rs = (360/DROT+DMIR);	
    int kn = k/(rs); 
    int rt = k%(rs); 
    
    bigw = nsize*(360/DROT+ DMIR);
		
	node_t cnode = dnodes[kn];
        int mid = bsize/2;
        	
        int rx = cnode.x-bsize/2;
        int ry = cnode.y-bsize/2;
        int rot = rt*DROT;
        

            cost_t c(node_t(rx,ry),0,rot,0,0);
            {
                //cout<<"hi1\n";
                c.cnode = dnodes[kn];
                //cout<<"hi2\n";
            	int kpos = (kn*rs)+rt;

                c.lpos = dem_lsizes[kpos].x;
                c.lsize = dem_lsizes[kpos].y;
                c.vpos = kpos*3;

                int cx = kpos; //colaesced memory access
             
                c.tnode = node_t(cx,0);
                c.skip  =false;


               for (int j=0; j<bsize; j++)
                    for (int i=0; i<bsize; i++)
                    {
                        //cout<<cx+i<<"/"<<cy+i<<": "<<cand(i,j)<<endl;
                        int ri = 0;
                        int rj = 0;
                        float candv=0;
			if (rt<360/DROT){
				 float ni=0;
                        	 float nj=0;
		                ni = (rx+mid) + ((i - mid)*cos_int(rot)) + ((j - mid) * sin_int(rot));
		                nj = (ry+mid) - ((i - mid)*sin_int(rot)) + ((j - mid) * cos_int(rot));
		                if (ni>=0 && nj>=0 && ni<src_w && nj<src_h)
		              		candv = cubicInterpol(src_ptr,src_w,src_h,ni,nj);    
                            }
                        else if (rot==360/DROT){
                        	ri = (rx+bsize-i-1);
                            	rj = (ry+j);
                            	if( ri>=0 && rj>=0 && ri<src_w && rj<src_h)
		                 		candv = src_ptr[ ri + rj*src_w  ];
                        }
                        else {
                            ri = (rx+i);
                            rj = (ry+bsize-j-1);
                             if( ri>=0 && rj>=0 && ri<src_w && rj<src_h)
                         		candv = src_ptr[ ri + rj*src_w  ];
                        
                        }                            
                                                 
                        
                        (bigtex)[cx] = candv;
                        cx += bigw; //colaesced memory access
                        if (candv<0.0001)   c.skip = true;
                    }
                    
                candidates[kpos] = (c);
            }
    }
}

void match_Feature_bef(Terrain& dest, Tree& usr_features, Tree& dem_features, vector<Image>& tar_pyr, vector<Image>&  src_pyr, int bsize)
{
    clock_t start_t, end_t, s_tmp, e_tmp;
    start_t = clock();
    
    match_time = 0;
    paste_time = 0;
    get_target = 0;

    Image target = usr_features.msource;

    Image src = dem_features.msource;

    node_list dem_nodes = dem_features.processNodes;



   int nsize = dem_nodes.size();
   int rs = 360/DROT;

    node_t* dnodes = new node_t [nsize];
    {
	    int count=0;
	    for (node_list::const_iterator it = dem_nodes.begin(); it != dem_nodes.end(); it++ ){
	    	dnodes[count] = *it;
	    	count++;
	    }
    }


    int bigwidth  = nsize*(rs+DMIR)*bsize*bsize;
    int bigheight = 1;
    bigw = nsize*(rs+DMIR);
    
    Image big(bigwidth,bigheight);
    float* bigtex = big.getPixels();
    float* src_ptr = src.getPixels();
    
    {
    end_t = clock();
    float elapsed = ((float)( end_t - start_t )) /CLOCKS_PER_SEC;
    cout<<"Get pyramid elapsed time: "<<elapsed<<" s.\n";
    start_t = clock();
    }

   float* variances = get_noise_stats_Feature(src_pyr, dem_nodes, bsize);
   
   {
    end_t = clock();
    float elapsed = ((float)( end_t - start_t )) /CLOCKS_PER_SEC;
    cout<<"Get noise stats elapsed time: "<<elapsed<<" s.\n";
    start_t = clock();
    }
    //float* variances = get_noise_stats_Feature(src_pyr, dem_nodes, bsize);

    //cost_t* candidates = new cost_t[nsize*rs];
    float* dem_vars = new float [nsize*(rs+DMIR)*3];
    node_t* dem_leafs = new node_t [nsize*(rs+DMIR)*10];
    node_t* dem_lsizes = new node_t [nsize*(rs+DMIR)];
    {
        int rt = 0;
        int rl = 0;
        int rv = 0;
        int cl=0;
        for (int kn=0; kn<nsize; kn++)
        {
            node_t cnode = dnodes[kn];
            for (int rot=0; rot<360; rot+=DROT)
            {
                vector<node_t> candleafs = getChildren(cnode,dem_features,bsize,rot);
                dem_lsizes[rl].x = cl;
                dem_lsizes[rl].y = candleafs.size();
                cl += dem_lsizes[rl].y; 
                rl++;
                for (unsigned int k=0; k<candleafs.size(); k++)
                    dem_leafs[rt++]=candleafs[k];
                for (int k=0; k<NLEVEL; k++)
                    dem_vars[rv++] = variances[kn*NLEVEL+k];
            }
            for (int m=0; m<DMIR; m++){
                int mir = m+1;
                vector<node_t> candleafs = getChildren(cnode,dem_features,bsize,0,mir);
                dem_lsizes[rl].x = cl;
                dem_lsizes[rl].y = candleafs.size();
                cl += dem_lsizes[rl].y;
                rl++;
                for (unsigned int k=0; k<candleafs.size(); k++)
                    dem_leafs[rt++]=candleafs[k];
                for (int k=0; k<NLEVEL; k++)
                    dem_vars[rv++] = variances[kn*NLEVEL+k];

              }
        }

    }
    
    {
    end_t = clock();
    float elapsed = ((float)( end_t - start_t )) /CLOCKS_PER_SEC;
    cout<<"Elapsed time: "<<elapsed<<" s.\n";
    start_t = clock();
    }

    float cand_time = 0;
    cout<<"Prepare candidates! "<< nsize<<"\n";
    start_t = clock();
    
    float* src_dev;		cudaMalloc((void**) &src_dev, sizeof(float)*src.getWidth()*src.getHeight());
    //float* dest_dev;		cudaMalloc((void**) &dest_dev, sizeof(float)*dest.getWidth()*dest.getHeight());
    //float* tar_dev;		cudaMalloc((void**) &tar_dev, sizeof(float)*target.getWidth()*target.getHeight());
    
    float* bigtex_dev;		cudaMalloc((void**) &bigtex_dev, sizeof(float)*bigwidth*bigheight);
    float* dem_vars_dev;	cudaMalloc((void**) &dem_vars_dev, sizeof(float)*nsize*(rs+DMIR)*3);
    
    thrust::device_vector<cost_t> candidates_dev(nsize*(rs+DMIR));
    //cost_t* candidates_dev;	cudaMalloc((void**) &candidates_dev, sizeof(cost_t)*nsize*rs);
    //cost_t* candidates_dev;
    //cudaMalloc((void**) &candidates_dev, sizeof)
    node_t* dnodes_dev;		cudaMalloc((void**) &dnodes_dev, sizeof(node_t)*nsize);
    
    node_t* dem_lsizes_dev;	cudaMalloc((void**) &dem_lsizes_dev, sizeof(node_t)*nsize*(rs+DMIR));
    node_t* dem_leafs_dev;	cudaMalloc((void**) &dem_leafs_dev, sizeof(node_t)*nsize*(rs+DMIR)*5);
    
    cudaMemcpy(src_dev, src.getPixels(), sizeof(float)*src.getWidth()*src.getHeight(), cudaMemcpyHostToDevice);
   // cudaMemcpy(dest_dev, dest.getPixels(), sizeof(float)*dest.getWidth()*dest.getHeight(), cudaMemcpyHostToDevice);
   // cudaMemcpy(tar_dev, target.getPixels(), sizeof(float)*target.getWidth()*target.getHeight(), cudaMemcpyHostToDevice);
    cudaMemcpy(bigtex_dev, big.getPixels(), sizeof(float)*big.getWidth()*big.getHeight(), cudaMemcpyHostToDevice);
    
    
    cudaMemcpy(dem_vars_dev, dem_vars, sizeof(float)*nsize*(rs+DMIR)*3, cudaMemcpyHostToDevice);
    //cudaMemcpy(candidates_dev, candidates, sizeof(cost_t)*nsize*rs, cudaMemcpyHostToDevice);
    cudaMemcpy(dnodes_dev, dnodes, sizeof(node_t)*nsize, cudaMemcpyHostToDevice);
    cudaMemcpy(dem_lsizes_dev, dem_lsizes, sizeof(node_t)*nsize*(rs+DMIR), cudaMemcpyHostToDevice);
    cudaMemcpy(dem_leafs_dev, dem_leafs, sizeof(node_t)*nsize*(rs+DMIR)*5, cudaMemcpyHostToDevice);
    
    int threadsPerBlock = 256;
    int blocksPerGrid =  ((nsize) / threadsPerBlock)+1;
	
    {	
    	dim3 dimBlock(32,rs+DMIR);
	dim3 dimGrid( (nsize / dimBlock.x)+1,1); //, A.height / dimBlock.y
	blocksPerGrid =  ((nsize*(rs+DMIR)) / threadsPerBlock)+1;
	
	//buildCandidates_Feature_Kernel<<<dimGrid, dimBlock>>> (thrust::raw_pointer_cast(&candidates_dev[0]), bigtex_dev, src_dev, dnodes_dev, dem_lsizes_dev,nsize, src.getWidth(),src.getHeight(),bsize);
	
	buildCandidates_Feature_Kernel<<<blocksPerGrid, threadsPerBlock>>> (thrust::raw_pointer_cast(&candidates_dev[0]), bigtex_dev, src_dev, dnodes_dev, dem_lsizes_dev,nsize, src.getWidth(),src.getHeight(),bsize);
    	
		
    	//cudaMemcpy(candidates, candidates_dev, sizeof(cost_t)*nsize*rs, cudaMemcpyDeviceToHost);
    	cudaMemcpy(big.getPixels(), bigtex_dev, sizeof(float)*bigwidth*bigheight, cudaMemcpyDeviceToHost);

    	end_t = clock();
    	float elapsed = ((float)( end_t - start_t )) /CLOCKS_PER_SEC;
    	cout<<"Elapsed time: "<<elapsed<<" s.\n";
    	cand_time += elapsed;
    	start_t = clock();
    }

    //big.savePGM("/tmp/bigtmp_Feature.pgm");	
    delete [] dem_lsizes;
    cudaFree(dem_lsizes_dev);

    node_list usr_nodes = usr_features.processNodes;
    cost_t* prev = new cost_t(node_t(-1,-1),0,0,0,0);
    
    float* usr_var_dev;    cudaMalloc((void**) &usr_var_dev,sizeof(float)*NLEVEL);
    node_t* uleafs_dev;    cudaMalloc((void**) &uleafs_dev,sizeof(node_t)*10);
    cost_t* prev_dev;      cudaMalloc((void**) &prev_dev,sizeof(cost_t));
    float* ucand_dev;	   cudaMalloc((void**) &ucand_dev,bsize*bsize*sizeof(float));
    float* utar_dev;       cudaMalloc((void**) &utar_dev,bsize*bsize*sizeof(float));
    
    cout<<"Start matching feature patches! "<<usr_nodes.size()<<" from "<<nsize*(rs+DMIR)<<" candidates\n";
    
    //return;
	
    int c_act = nsize*(rs+DMIR);

    int dx,dy;
    blocksPerGrid =  ((nsize*(rs+DMIR)) / threadsPerBlock)+1;
    for (node_list::const_iterator it = usr_nodes.begin(); it != usr_nodes.end(); it++ )
    {
        //dy=dx;
        s_tmp = clock();
        
        node_t cnode = *it;
        dx = cnode.x-bsize/2;
        dy = cnode.y-bsize/2;

        //cout<<"Leafs size: "<<usr_features.getcontrolpts(*it).size()<<endl;
        //target.get_crop(dx,dy,dx+(bsize-1),dy+(bsize-1),0,0).savePGM("/tmp/cand_usr.pgm");

        vector<node_t> uleafs = getChildren(cnode,usr_features,bsize);
        vector<float> usr_var =  noise_variances(tar_pyr,dx,dy,bsize);
        int lsize = imin(uleafs.size(),10);	

        Image ucand = dest.get_crop(dx,dy,dx+bsize-1,dy+bsize-1);
	Image utar = target.get_crop(dx,dy,dx+bsize-1,dy+bsize-1);
        
        
        cudaMemcpy(uleafs_dev, &uleafs[0], sizeof(node_t)*lsize, cudaMemcpyHostToDevice);
        cudaMemcpy(usr_var_dev, &usr_var[0], sizeof(float)*NLEVEL, cudaMemcpyHostToDevice);
	cudaMemcpy(ucand_dev, ucand.getPixels(), sizeof(float)*bsize*bsize, cudaMemcpyHostToDevice);
        cudaMemcpy(utar_dev, utar.getPixels(), sizeof(float)*bsize*bsize, cudaMemcpyHostToDevice);
	cudaMemcpy(prev_dev, prev, sizeof(cost_t), cudaMemcpyHostToDevice);

	//ComputeCosts_noFeature_Kernel<<<blocksPerGrid, threadsPerBlock>>>(ucand_dev, bigtex_dev,thrust::raw_pointer_cast(&candidates_dev[0]) , nsize*rs, dem_vars_dev, usr_var_dev, prev_dev, bsize);
	
	//cout<<"Look for help\n";
       
	ComputeCosts_Feature_Kernel<<<blocksPerGrid, threadsPerBlock>>>(ucand_dev, utar_dev, bigtex_dev,thrust::raw_pointer_cast(&candidates_dev[0]), nsize*(rs+DMIR),  dem_vars_dev, dem_leafs_dev, usr_var_dev, uleafs_dev, prev_dev, bsize, lsize);
	
	//cout<<"Got help\n";
	thrust::sort(candidates_dev.begin(),candidates_dev.end(),comp);	
	
	//cout<<"What now\n";

	//cudaMemcpy(candidates, candidates_dev, sizeof(cost_t)*nsize*rs, cudaMemcpyDeviceToHost);
        Image patch = findCand_dev(dest,bigtex,candidates_dev,nsize*(rs+DMIR), prev, bsize, dx,dy);
        //cout<<prev->org.x<<" "<<prev->org.y<<" "<<prev->rot<<" "<<prev->cost<<endl;
				if (prev->rot<360){
					Image patch2 = src.get_crop(prev->org.x,prev->org.y,prev->org.x+bsize-1,prev->org.y+bsize-1,prev->rot);
					//patch2.savePGM("./tmp/res_cand2.pgm");
					patch_merging(&dest, &patch2, dx, dy,1,bsize/10.);
				}
				else{
					Image patch2 = src.get_crop(prev->org.x,prev->org.y,prev->org.x+bsize-1,prev->org.y+bsize-1,0,(prev->rot-360)/DROT+1);
					//patch2.savePGM("./tmp/res_cand2.pgm");
					patch_merging(&dest, &patch2, dx, dy,1,bsize/10.);
				}
				//patch.savePGM("./tmp/res_cand1.pgm");
				//cin.get();
        
        //cudaMemcpy(prev_dev, prev, sizeof(cost_t), cudaMemcpyHostToDevice);
        
        e_tmp = clock();
        match_time+=mstimer(s_tmp,e_tmp);
        s_tmp = clock();
        
	//patch_merging(&dest, &patch, dx, dy,1,bsize/10.);
	
	e_tmp = clock();
        paste_time+=mstimer(s_tmp,e_tmp);
        //dest.savePGM("/tmp/res_tmp_cpu.pgm",dest.maxval);
        //dest.saveTerragen("/tmp/res_tmp_cpu.ter");

        //cin.get();

    }
    
    cudaFree(usr_var_dev);
    cudaFree(uleafs_dev);
    cudaFree(prev_dev);
     cudaFree(ucand_dev);
    cudaFree(utar_dev);
    delete prev;

    end_t = clock();
    float elapsed = ((float)( end_t - start_t )) /CLOCKS_PER_SEC;
    cout<<"Elapsed time: "<<elapsed<<" s.\n";

    //delete [] candidates;
    delete [] dnodes;
    delete [] dem_vars;
    delete [] variances;
    delete [] dem_leafs;
    
    //cudaFree(candidates_dev);
    cudaFree(dnodes_dev);
    cudaFree(dem_vars_dev);
    cudaFree(dem_leafs_dev);
    cudaFree(src_dev);
    //cudaFree(dest_dev);
    //cudaFree(tar_dev);
    cudaFree(bigtex_dev);
    
    cerr<<"\n\n*********** Feature matching GPU*******************\n";
     print_times();
    cerr<<" Candidates set: "<<cand_time<<"s\n";
    cerr<<" Number  of targets: "<<usr_nodes.size()<<"\n";
    cerr<<" Number  of cands: "<<nsize*(rs+DMIR)<<"\n";
    cerr<<"*********** End *******************\n\n";
}
