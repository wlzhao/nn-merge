#include "jmerge.h"
#include "iodelegator.h"
#include "kgraph.h"
#include <iostream>
#include <boost/timer/timer.hpp>
#include "vmath.h"


using namespace std;
using namespace kgraph;




Jmerge::Jmerge()
{

}



Jmerge::Jmerge(SpaceInterface<float> *func)
{
    setdisfunc(func);
}

void Jmerge::setdisfunc(SpaceInterface<float> *func)
{
    distfunc = func->get_dist_func();
}

Jmerge::~Jmerge()
{
    delete []GraphOne;

    GraphOne = NULL;

    delete []f_mat;

    f_mat = NULL;
}

unsigned* Jmerge::buildKGraph(float *mata, const unsigned K, const unsigned N, const unsigned D)
{
    cout<<"1. One Sub graph are prepared...................k = "<<K<<endl;
    InitialMemory(mata, K, N, D);
    preSize = row * 5 / 10 ;

    cout<<"Using NN-Descent build the first k-NN graph....."<<endl;
    f_mat = cutMat(mata, 0, preSize);

    KGraph *newGraph  = new KGraph(distfunc);

    GraphOne = newGraph->buildKGraph(f_mat, K0, preSize, col);

    this->nCmps += newGraph->nCmps;

    float scrate = this->nCmps/((this->row/2.0)*(this->row-1.0));
    ///cout<<"scan-rate: "<<scrate<<endl;

    cout<<"2. Joint Merge process.........................:"<<endl;
    this->nCmps = 0;
    boost::timer::cpu_timer timer;
    InitSGraph(0, GraphOne, 0, preSize);
    InitialGraph();

    nnDescent();
    merge2one(preSize);

    auto times = timer.elapsed();
    cout<<"Merge time .....................................times = "<<(times.wall / 1e9)<<endl;

    scrate = this->nCmps/((this->row/2.0)*(this->row-1.0));
    cout<<"Merge-Scanning rate:............................rate  = "<<scrate<<endl;


    ///delete newGraph;
    return 0;

}


unsigned Jmerge::InitialGraph()
{
    unsigned  rndid = 0;
    std::mt19937 rng(time(NULL));
    unsigned cutSize   = K0 * 3 / 6 ;/// K0 * 5 / 6

    unsigned RemainSize = this->K0 - cutSize;

    size_t tt = 0;

    ///cout<<row - preSize<<endl;
    ///#pragma omp parallel for default(shared) schedule(dynamic, 100) reduction(+:tt)
    for (unsigned i = 0; i < preSize; i ++)
    {
        unsigned *seeds = new unsigned[128];
        vector<unsigned> tp_seeds;
        egraph->cpnGraph[i].resize(cutSize);
        this->getRndSds(rng, 2*cutSize, seeds, row - preSize);
        for (unsigned j = 0; j < cutSize; j ++)
        {
            Neighbor &nb  = egraph->knnGraph[i].pool[j+RemainSize];
            Neighbor &rnb = egraph->cpnGraph[i][j];

            rnb.dist = nb.dist;
            rnb.flag = nb.flag;
            rnb.gd   = nb.gd;
            rnb.id   = nb.id;
        }
        for(unsigned m = 0; m < 2*cutSize; m ++)
        {
            bool db_flag = false;
            if(tp_seeds.size() >= cutSize)
                break;
            for(unsigned l = 0; l < egraph->knnGraph[i].pool.size(); l ++)
            {
                Neighbor &nb  = egraph->knnGraph[i].pool[l];
                rndid = seeds[m] +  preSize;  ///sample from the second graph.
                if(rndid == nb.id)
                {
                    db_flag = true;
                    break;
                }
            }
            if(!db_flag)
                tp_seeds.push_back(rndid);
        }
        for (unsigned j = 0; j < cutSize; j ++)
        {
            Neighbor &nb  = egraph->knnGraph[i].pool[j+RemainSize];
            nb.id    = tp_seeds[j];
            nb.gd    = 1;
            nb.flag  = true;
            nb.dist  = distfunc(mat, nb.id, mat, i, col);
            ++tt;
        }
        tp_seeds.clear();
        sort(egraph->knnGraph[i].pool.begin(), egraph->knnGraph[i].pool.end(), Neighbor::ComNeighbor);
    }


    nCmps += tt;

    size_t cc = 0;
    ///#pragma omp parallel for default(shared) schedule(dynamic, 100) reduction(+:cc)
    for (unsigned i = preSize; i < row; i ++)
    {
        unsigned *seeds = new unsigned[128];
        egraph->knnGraph[i].pool.resize(K0);
        egraph->knnGraph[i].L = K0;
        getRndSds(rng, K0+1, seeds, row);
        unsigned L = 0;
        for (unsigned j = 0; j < K0 + 1; j ++)
        {
            if(seeds[j] == i||L >= K0) continue;
            Neighbor &nb  = egraph->knnGraph[i].pool[L++];
            nb.id    = seeds[j];
            nb.gd    = 0;
            if(nb.id   >= preSize)
            {
                nb.gd = 1;
            }
            nb.flag  = true;
            nb.dist  = distfunc(mat, nb.id, mat, i, col);

            ++cc;
        }
        sort(egraph->knnGraph[i].pool.begin(), egraph->knnGraph[i].pool.end(), Neighbor::ComNeighbor);
        delete []seeds;
    }

    nCmps += cc;

    return 1;
}

unsigned Jmerge::localJoin()
{
    //unsigned begin_idx = 0, nid = 0, oid = 0, nnid = 0;
    ///unsigned Inserted = 0, Compared = 0;
    float dst = 0;
    size_t cc = 0;

    ///#pragma omp parallel for default(shared) schedule(dynamic, 100) reduction(+:cc)
    for(unsigned n = 0; n < row; n ++)
    {
        egraph->knnGraph[n].join([&](Ind i, Ind j) {
            if(i.id != j.id)
            {
                if(i.Gid != j.Gid||i.Gid == 1||j.Gid == 1)
                {
                    float dst = distfunc(mat, i.id, mat, j.id, col);
                    cc++;
                    egraph->knnGraph[i.id].parallel_try_insert(j.Gid, j.id, dst);
                    egraph->knnGraph[j.id].parallel_try_insert(i.Gid, i.id, dst);
                }
            }
        });
    }
    nCmps += cc;
    return 0;
 }

void Jmerge::test()
{
    const char *src = "../../../data/sift100k.txt";
    const char *dst = "../../../data/sift100k_Jmerge.txt";
    const char *grd = "../../../data/sift_learn_ground_truth_30.txt";
    unsigned d  = 0, n  = 0;
    unsigned k  = 20;
    float *refmat = IODelegator::load_refSet(src, d, n);

    L2space l2;
    ///L1space l1;
    Jmerge *myMerge = new Jmerge(&l2);
    myMerge->buildKGraph(refmat, k, n, d);
    myMerge->saveKGraph(dst);
    Evaluator::getKGraphRecall(dst, grd, 1);
    Evaluator::getKGraphRecall(dst, grd, 10);

    delete myMerge;

}

