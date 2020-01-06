#include "abstractmerge.h"
#include "iodelegator.h"
#include "vmath.h"
#include <iostream>
#include <algorithm>
///using namespace kgraph;
using namespace std;

///#include <boost/timer/timer.hpp>



AbstractMerge::AbstractMerge()
{

}

AbstractMerge::~AbstractMerge()
{
    clearGraph(egraph->knnGraph);
    clearGraph(egraph->cpnGraph);

}

unsigned AbstractMerge::InitialMemory(float *mata, const unsigned k, const unsigned n, const unsigned d)
{
    this->K0 = k;
    this->row = n;
    this->col = d;
    this->mat = mata;
    egraph = new MergeType();
    egraph->knnGraph.resize(n);
    egraph->cpnGraph.resize(n);
    return 0;
}

float* AbstractMerge::cutMat(float *mata, const unsigned b_idx, const unsigned l_idx)
{
    unsigned size_n = l_idx - b_idx;

    float *cmat = new float[size_n * col];

    for(unsigned i = b_idx, l = 0; i < l_idx && l < size_n; i ++, l++)
    {
        unsigned loc = i * col;
        unsigned reloc = l *col;
        for(unsigned j = 0; j < col; j ++)
        {
            cmat[reloc + j] = mata[loc + j];
        }
    }
    return cmat;
}

void AbstractMerge::InitSGraph(const unsigned b_idx, unsigned *graph, unsigned Gid , unsigned size_data)
{
    unsigned half = row / 2;

    for(unsigned i = 0; i < size_data; i ++)
    {
        unsigned abs_loc = b_idx + i;
        unsigned loc     = i * K0;
        egraph->knnGraph[abs_loc].pool.resize(K0);
        egraph->knnGraph[abs_loc].L = K0;
        for(unsigned j = 0; j < K0; j ++)
        {
            Neighbor &nb = egraph->knnGraph[abs_loc].pool[j];
            nb.id = graph[loc + j] + b_idx;
            nb.dist = distfunc(mat, abs_loc, mat, nb.id , col);

            nb.flag = true;
            nb.gd   = Gid;
        }
    }

}

unsigned AbstractMerge::nnDescent()
{
    ///unsigned hsize = row / 2;
    unsigned smplNum = 100;
    do
    {
        ///total++;
        ///#pragma omp parallel for
        for(unsigned irow = 0; irow < row; irow++)
        {
            auto &newlist  = egraph->knnGraph[irow].newlst;
            auto &oldlist  = egraph->knnGraph[irow].oldlst;
            unsigned nnSz = egraph->knnGraph[irow].pool.size();
            for(unsigned idim = 0; idim < nnSz; idim++)
            {
                Neighbor &nb = egraph->knnGraph[irow].pool[idim];
                auto &nhood_o = egraph->knnGraph[nb.id];
                float maxDst = egraph->knnGraph[nb.id].pool[nnSz - 1].dist;
                if(nb.flag)
                {
                    Ind tmp;
                    tmp.Gid = nb.gd;
                    tmp.id  = nb.id;
                    newlist.push_back(tmp);

                    if( maxDst < nb.dist )
                    {
                        LockGuard guard(nhood_o.lock);
                        Ind tmp;
                        tmp.Gid = 0;
                        if(irow >= preSize)
                        {
                            tmp.Gid = 1;
                        }
                        tmp.id  = irow;
                        nhood_o.rnewlst.push_back(tmp);
                    }
                    nb.flag = false;
                }
                else
                {
                    Ind tmp;
                    tmp.Gid = nb.gd;
                    tmp.id  = nb.id;
                    oldlist.push_back(tmp);
                    if(maxDst < nb.dist)
                    {
                        LockGuard guard(nhood_o.lock);
                        Ind tmp;
                        tmp.Gid = 0;
                        if(irow >= preSize)
                        {
                            tmp.Gid = 1;
                        }
                        tmp.id  = irow;
                        nhood_o.roldlst.push_back(tmp);
                    }
                }
            }///for(idim)
        }///for(irow)

        for( unsigned irow = 0; irow < row; irow ++ )
        {
            vector<Ind> &newlist  = egraph->knnGraph[irow].newlst;
            vector<Ind> &oldlist  = egraph->knnGraph[irow].oldlst;
            vector<Ind> &rnewlist = egraph->knnGraph[irow].rnewlst;
            vector<Ind> &roldlist = egraph->knnGraph[irow].roldlst;


            ///if(irow < hsize)
            {
                random_shuffle(rnewlist.begin(), rnewlist.end());
                if(rnewlist.size() > smplNum)
                {
                    rnewlist.resize(smplNum);
                }

                random_shuffle(roldlist.begin(), roldlist.end());

                if(roldlist.size() > smplNum)
                {
                    roldlist.resize(smplNum);
                }
            }


            newlist.insert(newlist.end(), rnewlist.begin(), rnewlist.end());
            oldlist.insert(oldlist.end(), roldlist.begin(), roldlist.end());
        }///for(irow)

        this->localJoin();

        ///--------------------------///
        accumulator_set<float, stats<tag::mean>> delta;
         for( unsigned irow = 0; irow < row; irow ++ ){
            delta(EvaluateDelta(egraph->knnGraph[irow].pool, egraph->knnGraph[irow].L));
        }

        ///cout<<mean(delta)<<endl;
        if(mean(delta)< default_delta)
        {
            break;
        }

        for(unsigned irow = 0; irow < row; irow++)
        {
            egraph->knnGraph[irow].newlst.clear();
            egraph->knnGraph[irow].oldlst.clear();
            egraph->knnGraph[irow].rnewlst.clear();
            egraph->knnGraph[irow].roldlst.clear();
        }
    }while(1);

    return 0;
}

float AbstractMerge::EvaluateDelta (mergeDataType::Neighbors &pool, unsigned K) {
    unsigned c = 0;
    unsigned N = K;
    if (pool.size() < N) N = pool.size();
    for (unsigned i = 0; i < N; ++i) {
        if (pool[i].flag) ++c;
    }
    return float(c) / K;  //the fraction of true in pool's flags
}

 void AbstractMerge::clearGraph(vector< vector<Neighbor> > &graph)
{
     unsigned i = 0, sz = 0;
     sz  = graph.size();

     for(i = 0; i < sz; i++)
     {
         vector<Neighbor> &nbs = graph[i];
         nbs.clear();
     }

     graph.clear();
}

void AbstractMerge::clearGraph(vector<NNList> &graph)
{
    unsigned sz = graph.size();
    for(unsigned i = 0; i < sz; i ++)
    {
        graph[i].pool.clear();
        graph[i].newlst.clear();
        graph[i].oldlst.clear();
        graph[i].rnewlst.clear();
        graph[i].roldlst.clear();
    }
    graph.clear();
}

unsigned AbstractMerge::merge2one(const unsigned sz)
{
    unsigned i, m = 0, n = 0, k, d1, d2;

    vector<Neighbor> tmpNbs;

    for (i = 0; i < sz; i ++)
    {
       d1 = egraph->knnGraph[i].pool.size();
       d2 = egraph->cpnGraph[i].size();
       tmpNbs.resize(d1);
       m = n = 0;
       for(k = 0; k < d1; k++)
       {
           Neighbor &m1 = egraph->knnGraph[i].pool[m];
           Neighbor &n1 = egraph->cpnGraph[i][n];
           Neighbor &k1 = tmpNbs[k];
           if(m1.dist <= n1.dist)
           {
               k1.dist = m1.dist;
               k1.flag = m1.flag;
               k1.id   = m1.id;
               k1.gd   = m1.gd;
               m++;
           }else{
               k1.dist = n1.dist;
               k1.flag = n1.flag;
               k1.id   = n1.id;
               k1.gd   = n1.gd;
               n++;
           }
       }
       for(k = 0; k < d1; k++)
       {
           Neighbor &m1 = egraph->knnGraph[i].pool[k];
           Neighbor &k1 = tmpNbs[k];
           m1.dist = k1.dist;
           m1.flag = k1.flag;
           m1.id   = k1.id;
           m1.gd   = k1.gd;
       }
    }
    tmpNbs.clear();
    clearGraph(egraph->cpnGraph);
    return 1;
 }
void AbstractMerge::saveKGraph(const char *src)
{
    ofstream *out = new ofstream(src, ios::out);
    (*out)<<row<<" "<<K0<<endl;
    for(unsigned i = 0; i < row; i ++)
    {
        (*out)<<i<<" "<<K0<<" ";
        for(unsigned j = 0; j < K0; j ++)
        {
            (*out)<<egraph->knnGraph[i].pool[j].id<<" ";
        }
        (*out)<<endl;
    }
}
