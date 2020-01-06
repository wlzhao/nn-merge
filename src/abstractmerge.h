#ifndef ABSMERGE_H
#define ABSMERGE_H
#include "abstractgraph.h"
#include "pvetex.h"
#include "space_lp.h"
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/moment.hpp>
using namespace boost;
using namespace boost::accumulators;

using namespace mergeDataType;


class AbstractMerge : public AbstractGraph
{
protected:
    MergeType *egraph;
    accumulator_set<float, stats<tag::mean>> delta;
    unsigned preSize;
public:
    AbstractMerge();
    ~AbstractMerge();
    virtual unsigned InitialMemory(float *mata, const unsigned k, const unsigned n, const unsigned d);
    virtual unsigned InitialGraph() = 0;
    virtual unsigned* buildKGraph(float *mata, const unsigned K, const unsigned N, const unsigned D) = 0;
    virtual void saveKGraph(const char *dst);
    virtual unsigned nnDescent();
    virtual unsigned localJoin() = 0;
    virtual void setdisfunc(SpaceInterface<float> *func) = 0;
    void clearGraph(vector< vector<Neighbor> > &graph);
    void clearGraph(vector<NNList> &graph);
    float EvaluateDelta (mergeDataType::Neighbors &pool, unsigned K);
    void InitSGraph(const unsigned b_idx, unsigned *graph, unsigned Gid, unsigned size_data);
    void mergeGraph(const char *src, const char *dst, const unsigned K);
    unsigned merge2one(const unsigned sz);
    float *cutMat(float *mata, const unsigned b_idx, const unsigned l_idx);
    ///static void pctest();
};



#endif // ABSMERGE_H
