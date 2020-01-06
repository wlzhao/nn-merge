#include "kgraph.h"
#include "iodelegator.h"
#include "vmath.h"

#include <boost/timer/timer.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/moment.hpp>
using namespace boost;
using namespace boost::accumulators;

using namespace kgraph;
using namespace std;

///const float KGraph::default_delta = 0.002;

KGraph::KGraph(SpaceInterface<float> *func)
{
    setdisfunc(func);
}

KGraph::KGraph(DISTFUNC<float> dfunc)
{
    distfunc = dfunc;
}

KGraph::~KGraph()
{
    delete []mat;
    mat = NULL;
}


unsigned KGraph::InitialMemory(float *mata, const unsigned k, const unsigned n, const unsigned d)
{
    this->mat = mata;
    this->nCmps = 0;
    this->row   = n;
    this->col = d;
    this->K0  = k;
    return 0;
}

float KGraph::EvaluateDelta (Neighbors  &pool, unsigned K) {
    unsigned c = 0;
    unsigned N = K;
    if (pool.size() < N) N = pool.size();
    for (unsigned i = 0; i < N; ++i) {
        if (pool[i].flag) ++c;
    }
    return float(c) / K;  //the fraction of true in pool's flags
}


unsigned KGraph::InitialGraph()
{
    unsigned i = 0, j = 0, L = 0;
    float dst = 0;

    unsigned *seeds  = new unsigned[512];

    knnGraph.resize(row);

    std::mt19937 rng(time(NULL));

    size_t cc = 0;
    ///#pragma omp parallel for default(shared) schedule(dynamic, 100) reduction(+:cc)
    for( i = 0; i < row; i++)
    {
        knnGraph[i].pool.resize(this->K0);
        knnGraph[i].L = this->K0;

        this->getRndSds(rng, this->K0+1, seeds, row);

        L = 0;
        for (j = 0; j < this->K0+1;  j++)
        {
            Neighborx &nb = knnGraph[i].pool[L];
            if(seeds[j] == i) continue;
            nb.id   = seeds[j];
            nb.pre_Owner = nb.id;  ///initial themselves.
            nb.dist = distfunc(mat, i, mat, nb.id, col);
            this->nCmps++;
            nb.flag = true;
            L ++;
            if(L >= this->K0)
                break;
        }
        sort(knnGraph[i].pool.begin(), knnGraph[i].pool.end(), Neighborx::ComNeighbor);
    }
    delete [] seeds;
    seeds = NULL;

    return 0;
}

unsigned KGraph::nndescent()
{
    unsigned smplNum = 100;
    float maxDst = 0;
    do
    {
        ///#pragma omp parallel for
        for(unsigned irow = 0; irow < row; irow++)
        {
            auto &newlist  = knnGraph[irow].newlst;
            auto &oldlist  = knnGraph[irow].oldlst;
            unsigned nnSz = knnGraph[irow].pool.size();
            for(unsigned idim = 0; idim < nnSz; idim++)
            {
                Neighborx &nb = knnGraph[irow].pool[idim];
                auto &nhood_o = knnGraph[nb.id];
                maxDst = knnGraph[nb.id].pool[nnSz - 1].dist;
                if(nb.flag)
                {
                    newlist.push_back(nb.id);

                    if( maxDst < nb.dist )
                    {
                        LockGuard guard(nhood_o.lock);
                        nhood_o.rnewlst.push_back(irow);
                    }
                    nb.flag = false;
                }
                else
                {
                    oldlist.push_back(nb.id);
                    if(maxDst < nb.dist)
                    {
                        LockGuard guard(nhood_o.lock);
                        nhood_o.roldlst.push_back(irow);
                    }
                }
            }///for(idim)
        }///for(irow)


        for( unsigned irow = 0; irow < row; irow ++ )
        {
            vector<unsigned> &newlist  = knnGraph[irow].newlst;
            vector<unsigned> &oldlist  = knnGraph[irow].oldlst;
            vector<unsigned> &rnewlist = knnGraph[irow].rnewlst;
            vector<unsigned> &roldlist = knnGraph[irow].roldlst;


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

        this->join();

        ///--------------------------///
        accumulator_set<float, stats<tag::mean>> delta;
        for( unsigned irow = 0; irow < row; irow ++ )
        {
            delta(EvaluateDelta(knnGraph[irow].pool, knnGraph[irow].L));
        }
        if(mean(delta) < AbstractGraph::default_delta)
        {
            break;
        }


        for(unsigned irow = 0; irow < row; irow++)
        {
            knnGraph[irow].newlst.clear();
            knnGraph[irow].oldlst.clear();
            knnGraph[irow].rnewlst.clear();
            knnGraph[irow].roldlst.clear();
        }

    }while(1);

    return 0;
}

unsigned KGraph::join()
{
    size_t cc = 0;
    ///#pragma omp parallel for default(shared) schedule(dynamic, 100) reduction(+:cc)
    for(unsigned n = 0; n < row; n ++)
    {
            knnGraph[n].join([&](unsigned i, unsigned j) {
            if( i != j  ) ///&& i.owner != j.owner
            {
                float dst = distfunc(mat, i, mat, j, col);
                cc++;
                knnGraph[i].parallel_try_insert(n, j, dst);
                knnGraph[j].parallel_try_insert(n, i, dst);
            }
        });
    }

    nCmps += cc;
    return 0;
}

void KGraph::setdisfunc(SpaceInterface<float> *func)
{
    distfunc = func->get_dist_func();
}

unsigned* KGraph::buildKGraph(float *mata , const unsigned k, const unsigned n, const unsigned d)
{
    boost::timer::cpu_timer timer;
    InitialMemory(mata, k, n, d);
    InitialGraph();
    nndescent();
    auto times = timer.elapsed();

    cout<<"Construction time ..............................times = "<<(times.wall / 1e9)<<endl;
    unsigned *knngraph = new unsigned[n*k];

    for(unsigned i = 0; i < n; i ++)
    {
        unsigned loc = i * k;
        for(unsigned j = 0; j < k; j ++)
        {
            knngraph[loc + j] = knnGraph[i].pool[j].id;
            ///cout<<knnGraph[i].pool[j].pre_Owner<<" ";
        }
        ///cout<<endl;
    }

    /**
    double scrate = this->nCmps/((this->row/2.0)*(this->row-1.0));
    cout<<this->nCmps<<endl;
    cout<<"scan-rate: "<<scrate<<endl;
    */
    return knngraph;
}

void KGraph::saveKGraph(const char *dst)
{
    ofstream *out = new ofstream(dst, ios::out);

    (*out)<<knnGraph.size()<<" "<<knnGraph[0].pool.size()<<endl;
    for(unsigned i = 0; i < knnGraph.size(); i ++)
    {
        (*out)<<i<<" "<< knnGraph[i].pool.size()<<" ";
        for(unsigned j = 0; j < knnGraph[i].pool.size(); j ++)
        {
            (*out)<<knnGraph[i].pool[j].id<<" ";
        }
        (*out)<<endl;
    }
}

void KGraph::test()
{
    const char *src = "../../data/sift100k.txt";
    const char *dst = "../../data/sift100k_KGraph.txt";
    const char *grd = "../../data/sift_learn_ground_truth_30.txt";
    unsigned k = 20;
    unsigned n = 100000;

    L2space l2;

    AbstractGraph *mygraph = new KGraph(&l2);

    unsigned d = 0;
    float *mat = IODelegator::load_refSet(src, d, n);

    mygraph->buildKGraph(mat, k, n,d);
    mygraph->saveKGraph(dst);

    Evaluator::getKGraphRecall(dst,grd, 1);
    Evaluator::getKGraphRecall(dst,grd, 10);

    delete mygraph;
}




