#ifndef DTTS_PPA_H
#define DTTS_PPA_H

#include "dtts_image.h"
using namespace Dtts;

typedef map<node_t,bool> mark_t;


/**
     * Tree class. Simple tree representing terrain ridges/valleys :-)
     *
     * @author Flora Tasse
     * @version 1.0.0
     */
class Tree
{

public:

    Tree();

    /**
     * get all nodes
     */
    node_list getNodes();

    /**
     * find neighbours
     */
    node_list& neighbours(const node_t& pnode);
    node_list neighbours_(const node_t& pnode);

    /**
     * bfs traversal
     */
    node_list bfs_full(const node_t& pnode);
    node_list bfs_iso(const node_t& pnode);

    node_list processNodes;
    vector<int> paths;

    int nedges;



    /**
     * restricted bfs traversal
     */
    node_t bfs_short(const node_t& parent, const node_t& child, node_list& visited, float rad);
    node_list bfs_short2(const node_t& parent, const node_t& child, node_list& visited, float rad);

    /**
     * Quickly get control point of a node
     */
    node_list getcontrolpts(const node_t& pnode);

    /**
    * Run PPA algorithm
    */
    bool ridge_type;
    void shaunPPA(Image& src, bool ridges=true, int srate = 2, int plength = 7, int portion = 20);

    void runPPA(bool is_ridge=true, int plength=5, int radius = 25, int groupmin=10, int branchmin=40, int samplerate=5, float ethres=0.0);

    void runPPA(Image& src, bool is_ridge=true, int plength=7, int radius = 25,int groupmin=10, int branchmin=40, int samplerate=5, float ethres=0.0)
    {
        clock_t tstart = clock();
        /*msource = src;
        cout<<"\nRun PPA on image ...\n";
        runPPA(is_ridge,plength,radius,groupmin,branchmin,samplerate,ethres);*/
        shaunPPA(src,is_ridge,samplerate, plength, branchmin);

        clock_t tstop = clock();
        //getComponents();

        if (roots.size()>0)
            root = roots[0];

        //renderFeatures(250).savePGM("/tmp/user_tmp.pgm");

        cerr<< "\n\n***************** MST Feature extraction terminated ********************\n";
        if (is_ridge)
            cerr<< " Ridge extraction from a terrain of size (" <<msource.width()<<","<<msource.height()<<")\n";
        else
            cerr<< " Valley extraction from a terrain of size (" <<msource.width()<<","<<msource.height()<<")\n";
        cerr<< " Number of features extracted: "<<getNodes().size()<<endl;
        /*cerr<< " Roots: "<<endl;
        for (int k=0; k<roots.size(); k++){
            cerr<<"  "<<roots[k].x<<"/"<<roots[k].y<<endl;
        }*/
        cerr<< " Elapsed time: "<<(float(difftime(tstop,tstart))/CLOCKS_PER_SEC)<<" s\n";
        cerr<< "***************** End ********************\n\n";
    }


    void runPPA_bis(Image& src, bool is_ridge=true, int plength=7, int radius = 25,int groupmin=10, int branchmin=40, int samplerate=5, float ethres=0.0)
    {
        clock_t tstart = clock();
        msource = src;
        cout<<"\nRun PPA on image ...\n";
        runPPA(is_ridge,plength,radius,groupmin,branchmin,samplerate,ethres);
        //shaunPPA(src,is_ridge,samplerate, plength, branchmin);

        clock_t tstop = clock();
        //getComponents();

        if (roots.size()>0)
            root = roots[0];

        //renderFeatures(250).savePGM("/tmp/user_tmp.pgm");

        cerr<< "\n\n***************** PPA Feature extraction terminated ********************\n";
        if (is_ridge)
            cerr<< " Ridge extraction from a terrain of size (" <<msource.width()<<","<<msource.height()<<")\n";
        else
            cerr<< " Valley extraction from a terrain of size (" <<msource.width()<<","<<msource.height()<<")\n";
        cerr<< " Number of features extracted: "<<getNodes().size()<<endl;
        /*cerr<< " Roots: "<<endl;
        for (int k=0; k<roots.size(); k++){
            cerr<<"  "<<roots[k].x<<"/"<<roots[k].y<<endl;
        }*/
        cerr<< " Elapsed time: "<<(float(difftime(tstop,tstart))/CLOCKS_PER_SEC)<<" s\n";
        cerr<< "***************** End ********************\n\n";
    }


    void Smoothing();



    /**
     * Load from file
     */
    void loadFile(const char* fname);
    void loadFile(Image& src, const char* fname, bool is_ridge=true){

        clock_t tstart = clock();

        msource = src;
        ridge_type = is_ridge;
        loadFile(fname);


    if (roots.size()>0)
            root = roots[0];

        clock_t tstop = clock();

        cerr<< "\n\n***************** Feature extraction terminated ********************\n";
        if (is_ridge)
            cerr<< " Ridge extraction from a terrain of size (" <<msource.width()<<","<<msource.height()<<")\n";
        else
            cerr<< " Valley extraction from a terrain of size (" <<msource.width()<<","<<msource.height()<<")\n";
        cerr<< " Number of features extracted: "<<getNodes().size()<<endl;
        /*cerr<< " Roots: "<<endl;
        for (int k=0; k<roots.size(); k++){
            cerr<<"  "<<roots[k].x<<"/"<<roots[k].y<<endl;
        }*/
        cerr<< " Elapsed time: "<<(float(difftime(tstop,tstart))/CLOCKS_PER_SEC)<<" s\n";
        cerr<< "***************** End ********************\n\n";
    }



    /**
     * Save to file
     */
    void saveFile(const char* fname);

    void saveFile_branches(const char* fname, unsigned int deg);

    Image renderFeatures(float color=100.);

    /**
    * Get control points for all nodes
    */
    void compute_control_pts(const float radius);

    // Source image
    Image msource;

    node_t root, prevNode;

    vector<node_t> components;

    map< node_t,vector<node_list> > followers;
    vector<point_t> getHprofile(node_t p);
    map_t isofeatures;

    vector<node_t> roots;
    void getRoots(int bsize);
    void getComponents();
    void branchpath(node_t p, int bsize);
    node_t branchpath(node_t p, node_t child, int bsize);

    Image renderFollowers(Image patch, int cx, int cy, int bsize, float color);


private:

    /**
     * test if node exists
     */
    bool is_node(node_t pnode);

    /**
     * test if an edge between two nodes exists
     */
    bool is_edge(const node_t& pnode1, const node_t& pnode2);

    /**
    * Insert an edge
    */
    void add_edge(const node_t& pnode1, const node_t& pnode2);

    /**
    * Remove an edge
    */
    void rem_edge(const node_t& pnode1, const node_t& pnode2);

    /**
     * Remove a node
     */
    void rem_node(const node_t& pnode);

    /**
     * Get control points of a feature
     */
    node_list control_pts(const node_t& pnode, const float radius);

    node_list branchlength(const node_t& parent, const node_t& child, unsigned int minimum);
    node_list branchlength2(const node_t& parent, const node_t& child, unsigned int minimum);

    node_list bfs(const node_t& pnode);

    void cleanBranches(unsigned int branchmin);
    void cleanBranches2(unsigned int branchmin);

    void cleanGroups(unsigned int groupmin);

    // Terrain features
    map_t mtree;

    //Control points
    map_t cpts;
};

#endif
