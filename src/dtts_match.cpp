#include "dtts_match.incl"

void match_Feature_cpu1(Terrain& dest, Tree& usr_features, Tree& dem_features, vector<Image>& tar_pyr, vector<Image>&  src_pyr,  int bsize)
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
    float* src_ptr = src.getPixels();
    int ncand = 1;

    vector<cost_t> candidates(nsize*(rs+DMIR));
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


    float cand_time=0;
    buildCandidates_Feature(&candidates[0], bigtex, src_ptr, dnodes, dem_lsizes,nsize, ncand, src.getWidth(),src.getHeight(),bsize);

    {
        e_tmp = clock();
        float elapsed = ((float)( e_tmp - s_tmp )) /CLOCKS_PER_SEC;
        cand_time += elapsed;
    }

    //big.savePGM("/tmp/bigtmp_Feature.pgm");
    delete [] dem_lsizes;

    node_list usr_nodes = usr_features.processNodes;
    cost_t prev (node_t(-1,-1),0,0,0,0);

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

        s_tmp = clock();
        ComputeCosts_Feature(ucand.getPixels(),target.getPixels(),bigtex,target.width(),ncand,&candidates[0],nsize*(rs+DMIR) ,&dem_vars[0],&dem_leafs[0],&usr_var[0],&uleafs[0],bsize, uleafs.size(), dx,dy);
        e_tmp = clock();
        match_time+=mstimer(s_tmp,e_tmp);

        s_tmp = clock();
        stable_sort(candidates.begin(),candidates.end(), comp);
        e_tmp = clock();
        sort_time+=mstimer(s_tmp,e_tmp);

        s_tmp = clock();
        Image patch = findCand(dest, bigtex, ncand, candidates,prev, bsize, dx,dy);
        e_tmp = clock();
        find_time+=mstimer(s_tmp,e_tmp);

        s_tmp = clock();
        patch_merging_cpu(&dest, &patch, dx, dy,1,bsize/10.);
        e_tmp = clock();
        paste_time+=mstimer(s_tmp,e_tmp);

    }

    delete [] dnodes;
    delete [] dem_vars;
    delete [] dem_leafs;

    cerr<<"\n\n*********** Feature matching CPU 1*******************\n";
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

void match_Feature_cpu2(Terrain& dest, Tree& usr_features, Tree& dem_features, vector<Image>& tar_pyr, vector<Image>&  src_pyr,  int bsize)
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

    float* src_ptr = src.getPixels();

    vector<cost_t> candidates(nsize*(rs+DMIR));
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


    float cand_time=0;
    buildCandidates_Feature(&candidates[0], src_ptr, dnodes, dem_lsizes,src.getWidth(),src.getHeight(),nsize,bsize);

    {
        e_tmp = clock();
        float elapsed = ((float)( e_tmp - s_tmp )) /CLOCKS_PER_SEC;
        cand_time += elapsed;
    }

    //big.savePGM("/tmp/bigtmp_Feature.pgm");
    delete [] dem_lsizes;

    node_list usr_nodes = usr_features.processNodes;
    cost_t prev (node_t(-1,-1),0,0,0,0);

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

        s_tmp = clock();
        ComputeCosts_Feature(ucand.getPixels(),target.getPixels(),src.getPixels(), target.width(), src.width(), src.height(),&candidates[0],nsize*(rs+DMIR) ,&dem_vars[0],&dem_leafs[0],&usr_var[0],&uleafs[0],bsize, uleafs.size(), dx,dy);
        e_tmp = clock();
        match_time+=mstimer(s_tmp,e_tmp);

        s_tmp = clock();
        stable_sort(candidates.begin(),candidates.end(), comp);
        e_tmp = clock();
        sort_time+=mstimer(s_tmp,e_tmp);

        s_tmp = clock();
        Image patch = findCand(dest,src,candidates, prev, bsize, dx,dy);
        e_tmp = clock();
        find_time+=mstimer(s_tmp,e_tmp);

        s_tmp = clock();
        patch_merging_cpu(&dest, &patch, dx, dy,1,bsize/10.);
        e_tmp = clock();
        paste_time+=mstimer(s_tmp,e_tmp);

    }

    delete [] dnodes;
    delete [] dem_vars;
    delete [] dem_leafs;

    cerr<<"\n\n*********** Feature matching CPU 2*******************\n";
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

void match_noFeature_cpu1(Terrain& dest, Image& src, Image& target, node_list dem_nodes, vector<Image>& tar_pyr, vector<Image>&  src_pyr,  int bsize, int osize)
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
    float* src_ptr = src.getPixels();
    int ncand = 1;

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


    float cand_time=0;
    //(cost_t* candidates, float* bigtex, float* src_ptr, node_t* dnodes, int nsize, int ncand, int src_w, int src_h, int bsize)
    buildCandidates_noFeature(&candidates[0], bigtex, src_ptr, dnodes, nsize, ncand, src.getWidth(),src.getHeight(),bsize);

    {
        e_tmp = clock();
        float elapsed = ((float)( e_tmp - s_tmp )) /CLOCKS_PER_SEC;
        cand_time += elapsed;
    }

    cost_t prev (node_t(-1,-1),0,0,0,0);

    int cnum = 0;
    //cout<<"Start matching feature patches! "<<usr_nodes.size()<<" from "<<nsize*(rs+DMIR)<<" candidates\n";

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
        //cout<<"hi0"<<endl;
        s_tmp = clock();
        ComputeCosts_noFeature(ucand.getPixels(), bigtex, ncand, &candidates[0], nsize*(rs+DMIR), dem_vars, &usr_var[0],bsize);
        e_tmp = clock();
        match_time+=mstimer(s_tmp,e_tmp);
        //cout<<"hi1"<<endl;

        s_tmp = clock();
        stable_sort(candidates.begin(),candidates.end(),comp);
        e_tmp = clock();
        sort_time+=mstimer(s_tmp,e_tmp);
        //cout<<"hi2"<<endl;

        s_tmp = clock();
        Image patch = findCand(dest, bigtex, ncand, candidates,prev, bsize, dx,dy);
        e_tmp = clock();
        find_time+=mstimer(s_tmp,e_tmp);
        //cout<<"hi3"<<endl;

        s_tmp = clock();
        patch_merging_cpu(&dest, &patch, dx, dy,1,bsize/10.);
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
                        ComputeCosts_noFeature(ucand.getPixels(), bigtex, ncand, &candidates[0], nsize*(rs+DMIR), dem_vars, &usr_var[0],bsize);
                        e_tmp = clock();
                        match_time+=mstimer(s_tmp,e_tmp);

                        s_tmp = clock();
                        stable_sort(candidates.begin(),candidates.end(),comp);
                        e_tmp = clock();
                        sort_time+=mstimer(s_tmp,e_tmp);

                        s_tmp = clock();
                        Image patch = findCand(dest, bigtex, ncand, candidates,prev, bsize, dx,dy);
                        e_tmp = clock();
                        find_time+=mstimer(s_tmp,e_tmp);

                        s_tmp = clock();
                        patch_merging_cpu(&dest, &patch, dx, dy,1,bsize/10.);
                        e_tmp = clock();
                        paste_time+=mstimer(s_tmp,e_tmp);

                        cnum++;
                    }

            if (finish)	break;
        }

    }



    delete [] dnodes;
    delete [] dem_vars;

    cerr<<"\n\n*********** Non feature matching CPU 1*******************\n";
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

void match_noFeature_cpu2(Terrain& dest, Image& src, Image& target, node_list dem_nodes, vector<Image>& tar_pyr, vector<Image>&  src_pyr,  int bsize, int osize)
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

    float* src_ptr = src.getPixels();

    vector<cost_t> candidates(nsize*(rs+DMIR));
    node_t* dnodes = new node_t [nsize];

    s_tmp = clock();

    float* dem_vars =  matchingPrepocessing(src_pyr, dem_nodes, dnodes, bsize);

    {
        e_tmp = clock();
        float elapsed = ((float)( e_tmp - s_tmp )) /CLOCKS_PER_SEC;
        cerr<<"Get noise stats elapsed time: "<<elapsed<<" s.\n";
        s_tmp = clock();
    }


    float cand_time=0;
    buildCandidates_noFeature(&candidates[0], src_ptr, dnodes, src.getWidth(), src.getHeight(), nsize, bsize);

    {
        e_tmp = clock();
        float elapsed = ((float)( e_tmp - s_tmp )) /CLOCKS_PER_SEC;
        cand_time += elapsed;
    }

    cost_t prev (node_t(-1,-1),0,0,0,0);
    int dx = -5*bsize, dy=-5*bsize;

    int cnum = nsize*(rs+DMIR);
    cout<<"Start matching non-feature patches! from "<<cnum<<"\n";
    list<cost_t> omega;
    float* conf = new float [dest.width()*dest.height()];
    for (int x=0; x<dest.width(); x++)	for (int y=0; y<dest.height(); y++)
            if (dest(x,y)>BG)
                conf[x+y*dest.width()] = 1;
            else
                conf[x+y*dest.width()] = 0;

    cnum = 0;
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
        ComputeCosts_noFeature(ucand.getPixels(),src.getPixels(), src.width(), src.height(),&candidates[0],nsize*(rs+DMIR) ,&dem_vars[0],&usr_var[0],bsize);
        e_tmp = clock();
        match_time+=mstimer(s_tmp,e_tmp);

        s_tmp = clock();
        stable_sort(candidates.begin(),candidates.end(), comp);
        e_tmp = clock();
        sort_time+=mstimer(s_tmp,e_tmp);

        s_tmp = clock();
        Image patch = findCand(dest,src,candidates, prev, bsize, dx,dy);
        e_tmp = clock();
        find_time+=mstimer(s_tmp,e_tmp);

        s_tmp = clock();
        patch_merging_cpu(&dest, &patch, dx, dy,1,bsize/10.);
        e_tmp = clock();
        paste_time+=mstimer(s_tmp,e_tmp);

        cnum++;
        //cout<<prev.org.x<<" "<<prev.org.y<< prev.rot<<endl;

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
                        ComputeCosts_noFeature(ucand.getPixels(),src.getPixels(), src.width(), src.height(),&candidates[0],nsize*(rs+DMIR) ,&dem_vars[0],&usr_var[0],bsize);
                        e_tmp = clock();
                        match_time+=mstimer(s_tmp,e_tmp);

                        s_tmp = clock();
                        stable_sort(candidates.begin(),candidates.end(), comp);
                        e_tmp = clock();
                        sort_time+=mstimer(s_tmp,e_tmp);

                        s_tmp = clock();
                        Image patch = findCand(dest,src,candidates, prev, bsize, dx,dy);
                        e_tmp = clock();
                        find_time+=mstimer(s_tmp,e_tmp);

                        s_tmp = clock();
                        patch_merging_cpu(&dest, &patch, dx, dy,1,bsize/10.);
                        e_tmp = clock();
                        paste_time+=mstimer(s_tmp,e_tmp);

                        cnum++;
                    }

            if (finish)	break;
        }

    }


    delete [] dnodes;
    delete [] dem_vars;

    cerr<<"\n\n*********** Non Feature matching CPU 2*******************\n";
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
    match_noFeature_cpu2(dest, src, target, dem_nodes, tar_pyr, src_pyr, bsize, osize);
}

void match_Feature(Terrain& dest, Tree& usr_features, Tree& dem_features, vector<Image>& tar_pyr, vector<Image>&  src_pyr, int bsize)
{
    match_Feature_cpu2(dest, usr_features, dem_features, tar_pyr, src_pyr, bsize);
}

