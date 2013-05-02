 # include <iostream>
 # include   <string.h>
 # include   <stdlib.h>

class Vertex_t
 {
    public:
       int x;
       int y;
       int label;

    public:
       Vertex_t( )   {  }
       ~Vertex_t( )  {  }

       void SetVertex_t(const int,const int,const int);

       Vertex_t(const int,const int,const int);

       //! set to value
        const Vertex_t &operator =(const Vertex_t &v){
            x= v.x; y= v.y;
            label = v.label;
            return *this;
        }

        //! test for equality
        bool operator==(const Vertex_t &v) const{	return (x==v.x && y==v.y);	}

        //! test for inequality
        bool operator!=(const Vertex_t &v)  const{	return (x!=v.x || y!=v.y);	}
 };

 class Edge_t
 {
    public:
       Vertex_t V1;
       Vertex_t V2;
       float weight;

    public:
       Edge_t( )   { }
       ~Edge_t( )  { }


       bool operator<(const Edge_t &other) const{
            return weight < other.weight;
       }

       void SetEdge_t(const Vertex_t,const Vertex_t,const float);

       Edge_t( const Vertex_t,const Vertex_t,const float ) ;

       //! set to value
        const Edge_t &operator =(const Edge_t &v){
            V1= v.V1; V2= v.V2;
            weight = v.weight;
            return *this;
        }

        //! test for equality
        bool operator==(const Edge_t &v) const{	return (V1==v.V1 && V2==v.V2);	}

        //! test for inequality
        bool operator!=(const Edge_t &v)  const{	return (V1!=v.V1 || V2!=v.V2);	}
 };

 void Vertex_t::SetVertex_t(const int _x,const int _y, const int _l)
 {
    x=_x;
    y=_y;
    label = _l;
 }

Vertex_t::Vertex_t(const int _x,const int _y, const int _l)
 {
    x=_x;
    y=_y;
    label = _l;
 }

 void Edge_t::SetEdge_t(const Vertex_t _V1,const Vertex_t _V2,const float _weight)
 {
    V1=_V1;
    V2=_V2;

    weight=_weight;
 }

 Edge_t::Edge_t(const Vertex_t _V1,const Vertex_t _V2,const float _weight)
 {
    V1=_V1;
    V2=_V2;

    weight=_weight;
 }


struct data{
int name;
int size;
struct data *home;
};
typedef struct data mydata;

class makeunionfind
{
public :
    mydata* S;
public:
    makeunionfind(int n)
    {
        S = new mydata [n];
        for(int i=0;i<n;i++)
        {
           S[i].name=i+1;
           S[i].size=0;

           S[i].home=&S[i];
        }

    }

    ~makeunionfind(){delete [] S; }

    void myunion(int a, int b)
    {
         int sizea,sizeb;
         sizea=S[a-1].home->size;
         sizeb=S[b-1].home->size;
        if(sizea>sizeb)
        {
           (S[b-1].home)->home=S[a-1].home;
           S[a-1].size=sizea+sizeb;

        }
        else
        {
            (S[a-1].home)->home=S[b-1].home;
            S[b-1].size=sizea+sizeb;

        }

    }
    int myfind(int a)
    {
        mydata *temp,*temp2,*stoppoint;
        int result;
        temp2=&S[a-1];
        temp=S[a-1].home;
        while(temp2!=temp)
        {
              temp2=temp;
              temp=(*temp2).home;
        }
        result=temp2->name;
         stoppoint=temp2;
           temp2=&S[a-1];
       //path compression
         do{
               temp=temp2->home;
               (*temp2).home=stoppoint;
               temp2=temp;
         }while(temp2!=stoppoint);
         return result;
    }
};

int compare(const void *a,const void *b )
{
	Edge_t *a1,*b1;
	a1=(Edge_t*)a;
	b1=(Edge_t*)b;
    if(a1->weight>b1->weight)
    return -1;
    else if (a1->weight<b1->weight)
    return 1;
    else
    return 0;

}

Edge_t *kruskal2(Edge_t *e, int n,int m, int *size)
{

    Edge_t *ans = new Edge_t [m];
    makeunionfind _list(n);

	int (*comp)(const void *a,const void *b );
	int k=0;
    comp=compare;
    qsort((void*)e,m,sizeof(e[0]),comp);

    for(int i=0;i<m;i++)
    {
        int s,f;
        s=(e[i].V1).label;
        f=(e[i].V2).label;

        if(_list.myfind(s)==_list.myfind(f))
        {
            continue;
        }
        else
        {
          _list.myunion(s,f);
          ans[k]=e[i];   //cout<<s<<" "<<f<<endl;
          k++;
        }

    }
    *size=k;
    return ans;

}
