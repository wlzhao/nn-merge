#include "lynndescent.h"

#include<omp.h>
#include "iodelegator.h"
#include "evaluator.h"
#include "timer.h"
///#include "vmath.h"
#include <boost/timer/timer.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/moment.hpp>

///#include "distopt.h"

#include <algorithm>
#include <iostream>
#include <cstring>
#include <cassert>
#include <string>
#include <ctime>

const unsigned int LyNNDescent::N0  = 64;
const unsigned int LyNNDescent::M_  = 16;
const float  LyNNDescent::default_delta = 0.002;


using namespace boost;
using namespace boost::accumulators;


LyNNDescent::LyNNDescent()
{
    this->col     = 0;
    this->row     = 0;
    this->K0      = 0;
    this->nCmps   = 0;
    this->S0      = 0;
    this->hgraph  = NULL;
    this->mat     = NULL;
}

LyNNDescent::LyNNDescent(SpaceInterface<float> *func)
{
    this->col     = 0;
    this->row     = 0;
    this->K0      = 0;
    this->nCmps   = 0;
    this->S0      = 0;
    this->hgraph  = NULL;
    this->mat     = NULL;
    distfunc = func->get_dist_func();

}

void LyNNDescent::initMemory(const char *srcFn)
{
    this->mat = IODelegator::load_refSet(srcFn, this->col, this->row);
    cout<<"Matrix size .............................. "<<this->row<<"x"<<this->col<<endl;
    this->total = this->row;
}


unsigned LyNNDescent::insert2Topk(const unsigned int id, const unsigned Gid,const unsigned int nid, const float dst)
{
    //kn is the k-nns
   //kn is the k-nns
    unsigned i, j, topk;
    vector<Neighbor> &crntLst = hgraph->knnGraph[id].pool;
    i = topk = crntLst.size();

    if(crntLst[topk-1].dist < dst)
        return 0;

    if(id == nid)
    return 0;

    while(i > 0)
    {
        j = i-1;
        if(crntLst[j].dist <= dst) break;
        i = j;
    }
    unsigned l = i;
    while(l > 0)
    {
        j = l - 1;
        if(crntLst[j].dist < dst) break;
        if(crntLst[j].id == nid) return 0;
        l = j;
    }

    j = topk - 1;

    while(j > i)
    {
       crntLst[j].id   = crntLst[j-1].id;
       crntLst[j].dist = crntLst[j-1].dist;
       crntLst[j].flag = crntLst[j-1].flag;
       crntLst[j].gd   = crntLst[j-1].gd;
        --j;
    }

    crntLst[j].id   = nid;
    crntLst[j].dist = dst;
    crntLst[j].flag = true;
    crntLst[j].gd   = Gid;

    return 1;
}


void LyNNDescent::getRndSds(std::mt19937 &rng, const unsigned int k, unsigned int *addr, const unsigned int N)
{
    if (N == k)
    {
        for (unsigned i = 0; i < k; ++i)
        {
            addr[i] = i;
        }
        return;
    }
    for (unsigned i = 0; i < k; ++i)
    {
        addr[i] = rng() % (N - k);
    }
    sort(addr, addr + k);
    for (unsigned i = 1; i < k; ++i)
    {
        if (addr[i] <= addr[i-1])
        {
            addr[i] = addr[i-1] + 1;
        }
    }
    unsigned off = rng() % N;
    for (unsigned i = 0; i < k; ++i)
    {
        addr[i] = (addr[i] + off) % N;
    }

}

void LyNNDescent::init()
{
    unsigned i = 0, j = 0, L = 0;

    unsigned *seeds  = new unsigned[512];
    hgraph = new GraphType();
    hgraph->absId     = new unsigned[this->row];
    hgraph->knnGraph.resize(LyNNDescent::N0);
    hgraph->cpnGraph.resize(LyNNDescent::N0);
    hgraph->inserted  = new unsigned char[this->row];
    hgraph->begin_idx = 0;
    hgraph->lst_idx   = LyNNDescent::N0;

    std::mt19937 rng(time(NULL));

    ///#pragma omp parallel for default(shared) schedule(dynamic, 100) reduction(+:cc)
    for( i = 0; i < LyNNDescent::N0; i++)
    {
        hgraph->knnGraph[i].pool.resize(this->K0);
        hgraph->knnGraph[i].L = this->K0;
        this->getRndSds(rng, this->K0+1, seeds, LyNNDescent::N0);
        L = 0;
        for (j = 0; j < this->K0+1;  j++)
        {
            Neighbor &nb = hgraph->knnGraph[i].pool[L];
            if(seeds[j] == i) continue;
            nb.id   = seeds[j];
            nb.dist = distfunc(mat, i, mat, nb.id, col);
            this->nCmps++;
            nb.flag = true;
            nb.gd   = 0;
            L ++;
            if(L >= this->K0)
                break;
        }
        sort(hgraph->knnGraph[i].pool.begin(), hgraph->knnGraph[i].pool.end(),Neighbor::ComNeighbor);
        hgraph->inserted[i] = true;
    }

    this->lastId = LyNNDescent::N0;
    this->total  = this->total - LyNNDescent::N0;

    delete [] seeds;
    seeds = NULL;
}

unsigned LyNNDescent::getGraphId(const unsigned Id)
{
    if(Id >= hgraph->begin_idx && Id < hgraph->lst_idx)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

unsigned LyNNDescent::nnDescent(unsigned char verse)
{
    unsigned threshold = 0, N = 0;
    unsigned smplNum = 100;
    unsigned long maxCmp = 0;
    N = hgraph->knnGraph.size();
    threshold = 0.001 * this->K0 * N;
    ///Ind tmp;
    float maxDst = 0;
    threshold = 10;
    do
    {
        ///#pragma omp parallel for
        for(unsigned irow = 0; irow < N; irow++)
        {
            auto &newlist  = hgraph->knnGraph[irow].newlst;
            auto &oldlist  = hgraph->knnGraph[irow].oldlst;
            unsigned nnSz = hgraph->knnGraph[irow].pool.size();
            for(unsigned idim = 0; idim < nnSz; idim++)
            {
                Neighbor &nb = hgraph->knnGraph[irow].pool[idim];
                auto &nhood_o = hgraph->knnGraph[nb.id];
                maxDst = hgraph->knnGraph[nb.id].pool[nnSz - 1].dist;
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
                        tmp.Gid = getGraphId(irow);
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
                        tmp.Gid = getGraphId(irow);
                        tmp.id  = irow;
                        nhood_o.roldlst.push_back(tmp);
                    }
                }
            }///for(idim)
        }///for(irow)


        for( unsigned irow = 0; irow < N; irow ++ )
        {
            vector<Ind> &newlist  = hgraph->knnGraph[irow].newlst;
            vector<Ind> &oldlist  = hgraph->knnGraph[irow].oldlst;
            vector<Ind> &rnewlist = hgraph->knnGraph[irow].rnewlst;
            vector<Ind> &roldlist = hgraph->knnGraph[irow].roldlst;


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

            newlist.insert(newlist.end(), rnewlist.begin(), rnewlist.end());
            oldlist.insert(oldlist.end(), roldlist.begin(), roldlist.end());
        }///for(irow)

        this->localJoin(maxCmp, verse);

        ///--------------------------///
        accumulator_set<float, stats<tag::mean>> delta;
         for( unsigned irow = 0; irow < N; irow ++ ){
            delta(EvaluateDelta(hgraph->knnGraph[irow].pool, hgraph->knnGraph[irow].L));
        }

        if(mean(delta)<LyNNDescent::default_delta)
        {
            break;
        }


        for(unsigned irow = 0; irow < N; irow++)
        {
            hgraph->knnGraph[irow].newlst.clear();
            hgraph->knnGraph[irow].oldlst.clear();
            hgraph->knnGraph[irow].rnewlst.clear();
            hgraph->knnGraph[irow].roldlst.clear();
        }
    }while(1);


    return 0;
}


unsigned LyNNDescent::localJoin(unsigned long &coutNN, unsigned char verse)
{
    //unsigned i = 0, j = 0, N = 0,l = 0, r = 0;
    //unsigned begin_idx, nid, oid, nnid;
    //unsigned Inserted = 0, Compared = 0;
    unsigned int N = 0;
    //float dst = 0;
    N      = hgraph->knnGraph.size();
    coutNN ++;
    size_t cc = 0;

    ///#pragma omp parallel for default(shared) schedule(dynamic, 100) reduction(+:cc)
    for(unsigned n = 0; n < N; n ++)
    {
        hgraph->knnGraph[n].join([&](Ind i, Ind j) {
            if(i.id != j.id)
            {
                if(i.Gid != j.Gid||i.Gid == 1||j.Gid == 1||verse)
                {
                    float dst = distfunc(mat, i.id, mat, j.id, col);
                    cc++;
                    hgraph->knnGraph[i.id].parallel_try_insert(j.Gid, j.id, dst);
                    hgraph->knnGraph[j.id].parallel_try_insert(i.Gid, i.id, dst);
                }
            }
        });
    }
    nCmps += cc;
    return 0;
 }

float LyNNDescent::EvaluateDelta (Neighbors const &pool, unsigned K) {
    unsigned c = 0;
    unsigned N = K;
    if (pool.size() < N) N = pool.size();
    for (unsigned i = 0; i < N; ++i) {
        if (pool[i].flag) ++c;
    }
    return float(c) / K;  //the fraction of true in pool's flags
}

 unsigned LyNNDescent::appdNwSmpl(const unsigned half, const unsigned bNum)
{
    unsigned i = 0, j = 0, sz = 0, hsize = 0, nk = 0, rndid = 0;
    std::mt19937 rng(time(NULL));
    unsigned *seeds = new unsigned[512];

    sz = hgraph->knnGraph.size();

    hgraph->cpnGraph.resize(sz);
    hgraph->lst_idx = this->lastId;

    vector<unsigned> tp_seeds;

    hsize = this->K0 - half;

    ///size_t cc = 0;
    ///#pragma omp parallel for default(shared) schedule(dynamic, 100) reduction(+:cc)
    for (i = 0; i < sz; i ++)
    {

        hgraph->cpnGraph[i].resize(hsize);
        this->getRndSds(rng, 2*half, seeds, bNum);
        for (j = 0; j < hsize; j ++)
        {
            Neighbor &nb  = hgraph->knnGraph[i].pool[j+half];
            Neighbor &rnb = hgraph->cpnGraph[i][j];
            rnb.dist = nb.dist;
            rnb.flag = nb.flag;
            rnb.gd   = nb.gd;
            rnb.id   = nb.id;
        }
        for(unsigned m = 0; m < 2*half; m ++)
        {
            bool db_flag = false;
            if(tp_seeds.size() >= hsize)
                break;
            for(unsigned l = 0; l < hgraph->knnGraph[i].pool.size(); l ++)
            {
                Neighbor &nb  = hgraph->knnGraph[i].pool[l];
                rndid = seeds[m] +  hgraph->lst_idx;
                if(rndid == nb.id)
                {
                    db_flag = true;
                    break;
                }
            }
            if(!db_flag)
                tp_seeds.push_back(rndid);
        }
        for (j = 0; j < hsize; j ++)
        {
            Neighbor &nb  = hgraph->knnGraph[i].pool[j+half];
            nb.id    = tp_seeds[j];
            nb.gd    = 1;
            nb.flag  = true;
            nb.dist = distfunc(mat, nb.id, mat, i, col);
            ++nCmps;
        }

        tp_seeds.clear();
        sort(hgraph->knnGraph[i].pool.begin(), hgraph->knnGraph[i].pool.end(), Neighbor::ComNeighbor);
    }

    nk = 0;
    hgraph->knnGraph.resize(sz + bNum);

    ///#pragma omp parallel for default(shared) schedule(dynamic, 100) reduction(+:cc)
    for(i = sz; i < hgraph->knnGraph.size(); i++)
    {
        nk = 0;
        hgraph->knnGraph[i].pool.resize(this->K0);
        hgraph->knnGraph[i].L = this->K0;
        this->getRndSds(rng, 2*this->K0, seeds, hgraph->knnGraph.size());

        for (j = 0; j < 2*this->K0; j++)
        {
            if(nk >= this->K0)
                break;
            Neighbor &nb = hgraph->knnGraph[i].pool[nk];
            if(seeds[j] == i) continue;
            nb.id   = seeds[j];
            nb.dist = distfunc(mat, i, mat, nb.id, col);
            nb.flag = true;
            nb.gd   = getGraphId(nb.id);
            nk ++;

            this->nCmps++;
        }

        sort(hgraph->knnGraph[i].pool.begin(), hgraph->knnGraph[i].pool.end(),Neighbor::ComNeighbor);
        hgraph->inserted[i] = true;
    }

    this->lastId = this->lastId  + bNum;
    this->total = this->total - bNum;

    delete [] seeds;
    return 1;
 }


unsigned LyNNDescent::merge2one()
{
    unsigned i, m = 0, n = 0, k = 0, sz = 0, d1 = 0, d2 = 0;
    sz = hgraph->cpnGraph.size();
    vector<Neighbor> tmpNbs;

    for (i = 0; i < sz; i ++)
    {
       d1 = hgraph->knnGraph[i].pool.size();
       d2 = hgraph->cpnGraph[i].size();
       tmpNbs.resize(d1);
       m = n = 0;
       for(k = 0; k < d1; k++)
       {
           Neighbor &m1 = hgraph->knnGraph[i].pool[m];
           Neighbor &n1 = hgraph->cpnGraph[i][n];
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
           Neighbor &m1 = hgraph->knnGraph[i].pool[k];
           Neighbor &k1 = tmpNbs[k];
           m1.dist = k1.dist;
           m1.flag = k1.flag;
           m1.id   = k1.id;
           m1.gd   = k1.gd;
       }
    }
    tmpNbs.clear();
    LyNNDescent::clearGraph(hgraph->cpnGraph);
    return 1;
 }

 unsigned LyNNDescent::dpg(ofstream *outStrm)
{
    unsigned i, j, sz, loc;

    sz  = hgraph->knnGraph.size();

    unsigned *knngraph = new unsigned[sz * this->K0];

    for(i = 0; i < sz; i ++)
    {
        loc = i * this->K0;
        for (j = 0; j < this->K0; j ++)
        {
            Neighbor &nb = hgraph->knnGraph[i].pool[j];
            knngraph[loc + j] =  nb.id;
        }
    }


    unsigned M_max = this->K0;
    unsigned M = this->K0/2;
    unsigned maxM = 0;

    if(sz == this->row) /// the base layer.
    {
        maxM = _M;
        ///maxM = 12;
    }
    else   ///not the base layer.
    {
        maxM = _M/2;
        ///maxM = 6;
    }

    this->DelaunayGraph(knngraph, sz, this->K0, M, maxM);

    this->saveDPGraph(outStrm);


    delete [] knngraph;

    return 0;
}

unsigned LyNNDescent::insert2topk(const unsigned id, const unsigned nid, const float dst, const unsigned maxSize)
{
    //kn is the k-nns
    unsigned i, j, topk;
    vector<Neighbor> &crntLst = this->graph[id];
    topk = crntLst.size();

    if(topk > 0)
    {
        i = topk;

        if(crntLst[topk-1].dist < dst && topk == maxSize )
            return 0;

        if(topk < maxSize)
        {
            if(crntLst[topk - 1].dist < dst)
            {
                crntLst.resize(i + 1);
                crntLst[topk].id   = nid;
                crntLst[topk].dist = dst;
                return 1;
            }
        }


        while (i > 0)
        {
            j = i - 1;
            if (crntLst[j].dist <= dst) break;
            i = j;
        }

        unsigned l = i; ///location///

        while(l > 0)
        {
            j = l - 1;
            if(crntLst[j].dist < dst) break;
            if(crntLst[j].id == nid) return 0;
            l = j;
        }


        if(topk < maxSize)
        {
            crntLst.resize(topk + 1);
            j = topk;
            while(j > i)
            {
               crntLst[j].id   = crntLst[j-1].id;
               crntLst[j].dist = crntLst[j-1].dist;
                --j;
            }
        }
        else
        {
            j = maxSize - 1;
            while(j > i)
            {
               crntLst[j].id   = crntLst[j-1].id;
               crntLst[j].dist = crntLst[j-1].dist;
                --j;
            }
        }

        crntLst[j].id   = nid;
        crntLst[j].dist = dst;
    }
    else
    {
        this->graph[id].resize(1);
        Neighbor &nb = this->graph[id][0];
        nb.dist = dst;
        nb.id   = nid;
    }
    return 1;
}

void LyNNDescent::saveDPGraph(ofstream *outStm)
{
    unsigned i , j;
    (*outStm)<< this->graph.size()<<endl;
    for(i = 0; i < this->graph.size(); i ++)
    {
        (*outStm) << i <<" "<<this->graph[i].size()<<" ";
        for(j = 0; j < this->graph[i].size(); j ++)
        {
            (*outStm)<<this->graph[i][j].id<<" ";
        }
        (*outStm)<<endl;
    }
}

void LyNNDescent::DelaunayGraph(unsigned *kNNGraph, const unsigned int N, const unsigned int K, const unsigned M_, const unsigned maxM)
{
    unsigned loc = 0, i = 0, j = 0, l = 0, n_b = 0;
    float dst = 0.0f, curdst = 0.0f;
    this->graph.resize(N);

    bool good = true;

    for(i = 0; i < N; i ++)
    {
        loc = i * K;
        for(j = 0; j < K; j ++)
        {
            if( this->graph[i].size() > M_ )
                break;

            n_b = kNNGraph[loc + j];

            dst = distfunc(this->mat, i, this->mat, n_b, this->col);
            good = true;

            for(l = 0; l < this->graph[i].size(); l ++)
            {
                curdst = distfunc(this->mat, graph[i][l].id, this->mat, n_b, this->col);
                if(curdst < dst)
                {
                    good = false;
                }

            }
            if(good)
            {
                insert2topk(i, n_b, dst, M_);
            }
        }
        ///cout<<"\r\r\r\r\t\t"<<i;
    }

    for(i = 0; i < N; i ++)
    {
        size_t NUM = this->graph[i].size();
        for(j = 0; j < NUM; j ++)
        {
            Neighbor &nn = this->graph[i][j];
            dst = nn.dist;

            if(this->graph[nn.id].size() < maxM)
            {
               insert2topk(nn.id, i, nn.dist, maxM);
            }

            else
            {
                bool  flag = true;
                for (size_t t = 0; t < this->graph[nn.id].size(); t ++)
                {
                    curdst = distfunc(this->mat, this->graph[nn.id][t].id, this->mat, i, this->col);
                    if(curdst < dst)
                    {
                        flag = false;
                        break;
                    }
                }
                if(flag)
                {
                    insert2topk(nn.id, i, nn.dist, maxM);
                }
            }

        }
        ///cout<<"\r\r\r\r\t\t"<<i;
    }
    return;
}



unsigned LyNNDescent::setflag()
{
     unsigned i = 0, j = 0, sz = 0, dim = 0;

     sz  = hgraph->knnGraph.size();
     for(i = 0; i < sz; i ++)
     {
         dim = hgraph->knnGraph[i].pool.size();
         for (j = 0; j < dim; j ++)
         {
             Neighbor &nb = hgraph->knnGraph[i].pool[j];
             nb.gd   = 0;
         }
     }
     return 0;
 }


 unsigned LyNNDescent::buildhgraph(const char *srcFn,  const char *dstFn, const char *grd, const unsigned topk, const unsigned M)
{
    string temp = dstFn;
    string SbDst = temp + "_baselayer.txt";
    string ShDst = temp + ".txt";
    ofstream *outStrm = new ofstream(ShDst, ios::out);
    const char *baseDst = SbDst.data();
    _M = M;
    //unsigned int i = 0, j = 0;
    //const unsigned Graph_num_ = 5,r = 20;
    unsigned char Tflag = 0, Fflag = 0;
    unsigned count_level = 0;
    long unsigned subtime = 0;
    this->K0 = topk;
    this->smplNum = LyNNDescent::N0;


    cout<<"Begin K-NN Graph construction ............ K="<<topk<<endl;
    boost::timer::cpu_timer timer;

    clock_t t1 = clock();
    this->initMemory(srcFn);
    this->init();

    Tflag = true;
    Fflag = false;
    this->Gid  = 1;
    this->nnDescent(Tflag);

    this->scrate = this->nCmps/((this->row/2.0)*(this->row-1.0));

    clock_t t2 = clock();
    this->dpg(outStrm);
    subtime += clock() - t2;

    count_level++;
    unsigned level[10] ={1,0,0,1,0,0,1,0,0,1};
    while(this->total != 0)
    {
        if(this->total > this->smplNum)
        {
            this->appdNwSmpl(this->K0/2, this->smplNum);
            this->nnDescent(Fflag);
            this->merge2one();

            if(level[count_level]&&count_level < 10)
            {
                clock_t t3 = clock();
                this->dpg(outStrm);
                subtime += clock() - t3;
            }
            count_level++;
            this->smplNum = hgraph->knnGraph.size();
            this->setflag();
            this->Gid++;
        }
        else
        {
            this->appdNwSmpl(this->K0/2, this->total);
            this->nnDescent(Fflag);
            this->merge2one();

            clock_t t4 = clock();
            this->dpg(outStrm);
            subtime += clock() - t4;

            count_level++;
            this->total = 0;
            this->Gid++;
        }

        this->scrate = this->nCmps/((this->row/2.0)*(this->row-1.0));

    }

    double t_time =  subtime/(CLOCKS_PER_SEC+0.0f);

    auto times = timer.elapsed();

    cout<<"Construction time ........................ times = "<<(times.wall / 1e9 - t_time)<<endl;

    this->scrate = this->nCmps/((this->row/2.0)*(this->row-1.0));

    cout<<"Total-comparions ......................... times = "<<this->nCmps<<endl;
    cout<<"Total-Scanning rate ...................... rate  = "<<this->scrate<<endl;


    this->saveKNNGraph(baseDst, topk);
    Evaluator::getKGraphRecall(baseDst, grd, 1);
    Evaluator::getKGraphRecall(baseDst, grd, 10);

    ///LyNNDescent::clearGraph(this->hgraph->knnGraph); ///need to be rewrite
    LyNNDescent::clearGraph(this->hgraph->cpnGraph);

    cout<<"save done!"<<endl;
    return 0;
}

void LyNNDescent::clearGraph(vector< vector<Neighbor> > &graph)
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

void LyNNDescent::saveKNNGraph(const char* dstFn, const unsigned int k0)
{
    unsigned i = 0, j = 0, k = k0;
    ofstream *outStrm = new ofstream(dstFn, ios::out);

    if(!outStrm->is_open())
    {
        cout<<"file: "<<dstFn<<" open failed!"<<endl;
        return ;
    }

    (*outStrm)<< this->hgraph->knnGraph.size()<<" "<<k0<<endl;

    for(i = 0; i < this->hgraph->knnGraph.size(); i ++)
    {
        k = k < this->hgraph->knnGraph[i].pool.size()?k:this->hgraph->knnGraph[i].pool.size();
        (*outStrm) << i <<" "<<k;
        for(j = 0; j < k; j++)
        {
            (*outStrm)<<" "<<this->hgraph->knnGraph[i].pool[j].id;
        }
        (*outStrm)<<endl;
    }

    outStrm->close();
}


LyNNDescent::~LyNNDescent()
{
    if(this->hgraph != NULL)
    {
       /// LyNNDescent::clearGraph(this->hgraph->knnGraph);
        LyNNDescent::clearGraph(this->hgraph->cpnGraph);
        delete [] this->hgraph->absId;
        delete [] this->hgraph->inserted;
        this->hgraph->absId    = NULL;
        this->hgraph->inserted = NULL;
        delete this->hgraph;
        this->hgraph = NULL;
    }

    if(this->mat != NULL)
    {
        delete [] this->mat;
        this->mat = NULL;
    }

    if(this->id2idx != NULL)
    {
        delete [] this->id2idx;
        this->id2idx = NULL;
    }
}


void LyNNDescent::test()
{

    const char *src1 = "../../data/sift100k.txt";
    const char *dst1 = "../../data/sift100k_HMGD_40";
    const char *grd1 = "../../data/sift_learn_ground_truth_30.txt";

    L2space l2;
    ///Ksquare ki;
    LyNNDescent *myknn = new LyNNDescent(&l2);

    myknn->buildhgraph(src1, dst1, grd1, 20, 20);

}
