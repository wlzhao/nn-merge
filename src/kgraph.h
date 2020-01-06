#ifndef KGRAPH_H
#define KGRAPH_H

#include "abstractgraph.h"
#include <boost/dynamic_bitset.hpp>
#include <unordered_set>
#include <vector>
#include "evaluator.h"
#include "space_lp.h"

#include<mutex>
#include "boost/smart_ptr/detail/spinlock.hpp"
//using namespace unordered_set;
using namespace std;


typedef boost::detail::spinlock Lock;
typedef std::lock_guard<Lock> LockGuard;



namespace kgraph{

        struct Neighborx {
            unsigned id;
            int pre_Owner;
            float dist;
            bool flag;  // whether this entry is a newly found one
            Neighborx () {}
            Neighborx (int pre, unsigned i, float d, bool f = true):pre_Owner(pre),id(i), dist(d), flag(f) {
            }
            static int ComNeighbor(Neighborx A, Neighborx B)
            {
                return A.dist < B.dist;
            }
        };

        typedef std::vector<Neighborx> Neighbors;



        static  unsigned UpdateKnnList (Neighborx *addr, unsigned K, Neighborx nn) {
        // find the location to insert

        unsigned i =0, topk =0, j = 0;

        i = topk = K;

        if(addr[topk-1].dist <= nn.dist)
            return 0;


        while(i > 0)
        {
            j = i-1;
            if(addr[j].dist <= nn.dist) break;
            i = j;
        }
        unsigned l = i;
        while(l > 0)
        {
            j = l - 1;
            if(addr[j].dist < nn.dist) break;
            if(addr[j].id == nn.id) return 0;
            l = j;
        }

        j = topk - 1;

        while(j > i)
        {
           addr[j].id   = addr[j-1].id;
           addr[j].dist = addr[j-1].dist;
           addr[j].flag = addr[j-1].flag;
           addr[j].pre_Owner = addr[j - 1].pre_Owner;
            --j;
        }

        addr[j].id   = nn.id;
        addr[j].dist = nn.dist;
        addr[j].flag = nn.flag;
        addr[j].pre_Owner = nn.pre_Owner;

        return 1;
    }


    struct NNList {
            unsigned L;
            Lock lock;
            vector<unsigned> newlst;
            vector<unsigned> oldlst;
            vector<unsigned> rnewlst;
            vector<unsigned> roldlst;
            Neighbors pool;
            float maxDst;

            unsigned parallel_try_insert (int pre , unsigned id, float dist) {

                    ///if (dist > radius) return pool.size();
                    LockGuard guard(lock);
                    //a neighbor of a neighbor is also likely to be a neighbor.
                    // where the location of the "id" to insert.
                    UpdateKnnList(&pool[0], L, Neighborx(pre ,id, dist, true));
                    //if the location is smaller than neighborhood size L it.
                    //
                    return 0;
                }

              template <typename C>
                void join (C callback) const {
                    for (unsigned i: newlst) {
                        for (unsigned j: newlst) {
                            if (i  < j ) {
                                callback(i, j);
                            }
                        }
                        for (unsigned j: oldlst) {
                            callback(i, j);
                        }
                    }
                }
            };


    class KGraph: public AbstractGraph{

    private:
        ///static const float default_delta;
    protected:

        vector<int> visit_list;
        std::vector<NNList> knnGraph;

    public:
            KGraph():AbstractGraph(){};
            KGraph(SpaceInterface<float> *func);
            KGraph(DISTFUNC<float> dfunc);
            ///KGraph(unsigned k, unsigned s): K0(k), row(0), col(0), SizeN(s), nCmps(0), mat(NULL){};
            ~KGraph();
            unsigned InitialMemory(float *src, const unsigned k, const unsigned n, const unsigned d);
            float    EvaluateDelta (Neighbors  &pool, unsigned K);
            unsigned InitialGraph();
            unsigned nndescent();
            unsigned join();
            static void test();
            virtual unsigned* buildKGraph(float *mata, const unsigned K, const unsigned N, const unsigned d);
            virtual void saveKGraph(const char *dst);
            virtual void setdisfunc(SpaceInterface<float> *func);
    };
}





#endif // KGRAPH_H
