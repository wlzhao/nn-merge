#include "smerge.h"
#include "iodelegator.h"
#include "kgraph.h"
#include <iostream>
#include <boost/timer/timer.hpp>
#include "vmath.h"


using namespace std;
using namespace kgraph;




Smerge::Smerge()
{

}

Smerge::Smerge(SpaceInterface<float> *func)
{
    ///AbstractGraph
    setdisfunc(func);
}

Smerge::~Smerge()
{
    delete []GraphOne;
    delete []GraphTwo;

    GraphOne = NULL;
    GraphTwo = NULL;

    delete []f_mat;
    delete []b_mat;

    f_mat = NULL;
    b_mat = NULL;
}

void Smerge::setdisfunc(SpaceInterface<float> *func)
{
    distfunc = func->get_dist_func();
}

unsigned* Smerge::buildKGraph(float *mata, const unsigned K, const unsigned N, const unsigned D)
{
    ///boost::timer::cpu_timer timer;

    cout<<"1. Two Sub graphs are prepared..................k = "<<K<<endl;
    InitialMemory(mata, K, N, D);

    preSize = row * 5 / 10; ///row * 2 / 8

    f_mat = cutMat(mata, 0, preSize);
    b_mat = cutMat(mata, preSize, row);

    cout<<"Using NN-Descent build the first k-NN graph....."<<endl;
    KGraph *newGraph  = new KGraph(distfunc);
    GraphOne = newGraph->buildKGraph(f_mat, K0, preSize, col);

    this->nCmps += newGraph->nCmps;

    cout<<"Using NN-Descent build the second k-NN graph...."<<endl;
    KGraph *newGraph1  = new KGraph(distfunc);
    GraphTwo = newGraph1->buildKGraph(b_mat, K0, row - preSize, col);

    this->nCmps += newGraph1->nCmps;

    float scrate = this->nCmps/((this->row/2.0)*(this->row-1.0));
    ///cout<<"sub-rate: "<<scrate<<endl;

    this->nCmps = 0;
    cout<<"2. Symmetric Merge process.....................:"<<endl;
    boost::timer::cpu_timer timer;
    InitSGraph(0, GraphOne, 0, preSize);
    InitSGraph(preSize, GraphTwo, 1, row - preSize);
    InitialGraph();
    nnDescent();
    merge2one(row);
    auto times = timer.elapsed();
    cout<<"Merge time .....................................times = "<<(times.wall / 1e9)<<endl;



    scrate = this->nCmps/((this->row/2.0)*(this->row-1.0));
    cout<<"Merge-Scanning rate:............................rate  = "<<scrate<<endl;

    return 0;

}


unsigned Smerge::InitialGraph()
{
    unsigned  rndid = 0;
    std::mt19937 rng(time(NULL));

    unsigned cutSize   = K0 * 1 / 2;///K0 * 1 / 6
    unsigned RemainSize = this->K0 - cutSize;

    size_t tt = 0;

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
        delete [] seeds;
        sort(egraph->knnGraph[i].pool.begin(), egraph->knnGraph[i].pool.end(), Neighbor::ComNeighbor);
    }
    nCmps += tt;



    size_t cc = 0;
    ///#pragma omp parallel for default(shared) schedule(dynamic, 100) reduction(+:cc)
    for (unsigned i = preSize; i < row; i ++)
    {
        unsigned *seeds = new unsigned[128];
        vector<unsigned> tp_seeds;
        egraph->cpnGraph[i].resize(cutSize);
        this->getRndSds(rng, 2*cutSize, seeds, preSize);
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
                rndid = seeds[m]; ///sample from the first graph
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
            nb.gd    = 0;
            nb.flag  = true;
            nb.dist  = distfunc(mat, nb.id, mat, i, col);
            ++cc;
        }
        tp_seeds.clear();
        delete []seeds;
        sort(egraph->knnGraph[i].pool.begin(), egraph->knnGraph[i].pool.end(), Neighbor::ComNeighbor);
    }

    nCmps += cc;

    return 1;
}


unsigned Smerge::localJoin()
{
    ///unsigned i = 0, j = 0, N = 0,l = 0, r = 0;
    ///unsigned begin_idx, nid, oid, nnid;
    ///unsigned Inserted = 0, Compared = 0;
    ///float dst = 0;
    size_t cc = 0;

    ///#pragma omp parallel for default(shared) schedule(dynamic, 100) reduction(+:cc)
    for(unsigned n = 0; n < row; n ++)
    {
        egraph->knnGraph[n].join([&](Ind si, Ind sj) {
            if(si.id != sj.id && si.Gid != sj.Gid)
            {
                float dst = distfunc(mat, si.id, mat, sj.id, col);
                cc++;
                egraph->knnGraph[si.id].parallel_try_insert(sj.Gid, sj.id, dst);
                egraph->knnGraph[sj.id].parallel_try_insert(si.Gid, si.id, dst);
            }
        });
    }
    nCmps += cc;
    return 0;
}


void Smerge::test()
{

    const char *src = "../../../data/sift100k.txt";
    const char *dst = "../../../data/sift100k_Smerge.txt";
    const char *grd = "../../../data/sift_learn_ground_truth_30.txt";
    unsigned d  = 0, n  = 0;
    unsigned k  = 20;


    float *refmat = IODelegator::load_refSet(src, d, n);

    L2space l2;
    ///L1space l1;

    AbstractGraph *myMerge = new Smerge(&l2);
    myMerge->buildKGraph(refmat, k, n, d);
    myMerge->saveKGraph(dst);
    Evaluator::getKGraphRecall(dst, grd, 1);
    Evaluator::getKGraphRecall(dst, grd, 10);

    delete myMerge;

}

