#ifndef LYNNDESCENT_H
#define LYNNDESCENT_H

#include "pvetex.h"
#include<iostream>
///#include "rtpair.h"

#include <boost/dynamic_bitset.hpp>
#include <unordered_set>
#include <vector>
#include "space_lp.h"

using namespace mergeDataType;





class LyNNDescent
{
    public:
        LyNNDescent();
        LyNNDescent(SpaceInterface<float> *func);
        virtual ~LyNNDescent();

    private:
    static const unsigned int N0;
    static const unsigned int M_;
    static const float default_delta;

protected:
    unsigned K0, S0;
    unsigned row, col, sizeN;
    unsigned long nCmps;
    unsigned smplNum;
    unsigned total;
    unsigned lastId, Gid;
    float *mat, scrate;
    unsigned char *g_flag;
    unordered_set<unsigned> vlst;
    DISTFUNC<float> distfunc;
    double total_time = 0;
    unsigned avgVisit, maxlevel_, _M;
    GraphType *hgraph;
    ///vector<NNList> nhood;

    ///for dpg
    unsigned *augKnnGraph;
    unsigned *knnLocs, *id2idx;
    unsigned idx_size;

    vector<Neighbors> graph;


public:

    void initMemory(const char *src);
    void init();

    unsigned nnDescent(unsigned char verse);
    unsigned getGraphId(const unsigned Id);
    unsigned dpg(ofstream *outStrm);
    void DelaunayGraph(unsigned *kNNGraph, const unsigned int N, const unsigned int K, const unsigned M_, const unsigned maxM);

    float EvaluateDelta (Neighbors const &pool, unsigned K);
    unsigned localJoin(unsigned long &coutNN, unsigned char verse);
    unsigned appdNwSmpl(const unsigned half, const unsigned bNum);
    void saveDPGraph(ofstream *outStm);
    void     getRndSds(std::mt19937 &rng, const unsigned int k, unsigned int *addr, const unsigned int N);
    unsigned insert2Topk(const unsigned int id, const unsigned Gid,const unsigned int nid, const float dst);
    unsigned insert2topk(const unsigned id, const unsigned nid, const float dst, const unsigned maxSize);
    ///unsigned insert2topk(const unsigned id, const unsigned nid, const unsigned dst, const unsigned maxSize);
    unsigned merge2one();
    unsigned setflag();
    unsigned buildhgraph(const char *srcFn,  const char *dstFn, const char *grd,const unsigned topk, const unsigned M);///entrance///
    static void  clearGraph(vector< vector<Neighbor> > &graph);

    static void test();

    ///---------------IO-----------------///
    void saveKNNGraph(const char* dstFn, const unsigned int k0);
    void saveKNNGraph(ofstream *outStrm, const unsigned int k0);

    ///void printTops(unsigned int i, RTPair *pairs, ostream *outStrm, int topk);

};

#endif // LYNNDESCENT_H
