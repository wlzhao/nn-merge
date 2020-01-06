#ifndef LYSEARCH_H
#define LYSEARCH_H

#include "rtpair.h"
#include "space_lp.h"

#include <boost/dynamic_bitset.hpp>
#include <unordered_set>
#include <vector>

using namespace std;


struct LyGraph
{
    unsigned *augKnnGraph;
    unsigned *knnLocs;
    unsigned ilayer;
    unsigned szGraph;
    unsigned upSz;
};

class LySearch
{
public:
    LySearch();
    LySearch(const char *conf);
    LySearch(const char *Rsrc, const char *Gsrc, SpaceInterface<float> *func);
    LySearch(SpaceInterface<float> *func);
    virtual ~LySearch();
    struct CompareBySecond {
    constexpr bool operator()(std::pair<int, float> const &a,
                              std::pair<int, float> const &b) const noexcept {
            return a.second < b.second;
        }
    };
    struct CompareByFirst {
            constexpr bool operator()(std::pair<float, unsigned> const &a,
                                      std::pair<float, unsigned> const &b) const noexcept {
                return a.first < b.first;
            }
        };



private:
    static const unsigned int Sd, L0;
    static const unsigned int maxlevel;

protected:
    unsigned K0, S0;
    unsigned row, col, sizeN, nQry, qDim;
    unsigned long nCmps;
    unsigned smplNum, nLayer;
    DISTFUNC<float> distfunc;
    unsigned short int *flags;
    float *mat, scrate, *queries, *crntQry;
    double total_time = 0, kk_time = 0;
    ///for dpg
    unsigned *crntAugKnnGraph;
    unsigned *crntKnnLocs;
    vector<unsigned int> pts;
    LyGraph *layers;

    int  initconfig(const char *refer, const char *Graph, const char * query);
    int  knnSearch(const char *qryFn, const char *dstFn, const char *grdFn,const char *record);
    int  singleknnSearch(const char *qryFn, const char *dstFn, const char *grdFn,const char *record);
    float nndFlatSrch(unsigned maxlayer, vector<RTPair> &results, const unsigned topk);

    int  nndWDSearch(vector<RTPair> &topLst, const unsigned int topk, const unsigned int S0);
    float nndHNSWSrch(unsigned cur_obj, float cur_dst, unsigned maxlayer, vector<RTPair> &results, const unsigned topk);
    int  loadRefSet(const char *srcFn);
    int  loadGraphs(const char *srcFn);
    int  loadsingleGraph(const char *srcFn);
    unsigned insert2List(RTPair *addr, unsigned K, RTPair &nwnn);
    unsigned insert2List(RTPair *addr, unsigned K, unsigned id, float dst, const unsigned Sz);
    void releaseGraphs(LyGraph *myLayers, const unsigned int nLayer);
    bool merge2one(vector<RTPair> &major, vector<RTPair> &minor, const unsigned topk);
    void getRndSds(std::mt19937 &rng, const unsigned int k, unsigned int *addr, const unsigned int N);
    void printTops(unsigned int i, vector<RTPair> &pairs, ostream *outStrm, int topk);
    void saveBasedlayer(const char *dst);




public:
    static void test();
    static void flatsearch();
    static void savebase();

};

#endif
