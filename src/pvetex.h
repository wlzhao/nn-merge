#ifndef PVETEX_H
#define PVETEX_H
#include <vector>

#include<mutex>
#include "boost/smart_ptr/detail/spinlock.hpp"
//using namespace unordered_set;
using namespace std;


typedef boost::detail::spinlock Lock;
typedef std::lock_guard<Lock> LockGuard;


namespace mergeDataType{


struct Neighbor {
    unsigned id;
    float dist;
    unsigned gd;
    bool flag;  // whether this entry is a newly found one
    Neighbor () {}
    Neighbor (unsigned g,unsigned i, float d, bool f = true):gd(g), id(i), dist(d), flag(f) {
    }
    static int ComNeighbor(Neighbor A, Neighbor B)
    {
        return A.dist < B.dist;
    }
};

typedef std::vector<Neighbor> Neighbors;


struct Ind {
        unsigned id;
        unsigned Gid;
};

static  unsigned UpdateKnnList (Neighbor *addr, unsigned K, Neighbor nn) {
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
           addr[j].gd   = addr[j-1].gd;
            --j;
        }

        addr[j].id   = nn.id;
        addr[j].dist = nn.dist;
        addr[j].flag = nn.flag;
        addr[j].gd   = nn.gd;

        return 1;

    }


struct NNList {
        unsigned L;
        Lock lock;
        vector<Ind> newlst;
        vector<Ind> oldlst;
        vector<Ind> rnewlst;
        vector<Ind> roldlst;
        Neighbors pool;
        float maxDst;


        unsigned parallel_try_insert (unsigned Gid, unsigned id, float dist) {

                ///if (dist > radius) return pool.size();
                LockGuard guard(lock);
                //a neighbor of a neighbor is also likely to be a neighbor.
                // where the location of the "id" to insert.
                unsigned l = UpdateKnnList(&pool[0], L, Neighbor(Gid, id, dist, true));
                //if the location is smaller than neighborhood size L it.
                //
                return 0;
            }

         template <typename C>
            void join (C callback)  const{
                unsigned nsize = newlst.size();
                unsigned osize = oldlst.size();

                for (unsigned i = 0; i < nsize; i++) {
                    for (unsigned j = i + 1;  j < nsize; j++) {
                        callback(newlst[i], newlst[j]);
                    }
                    for (unsigned j = 0;  j < osize; j++) {
                        callback(newlst[i], oldlst[j]);
                    }
                }

            }
};

    struct GraphType{
    public:
        unsigned *absId;
        unsigned char *inserted;
        std::vector<NNList> knnGraph;
        std::vector<Neighbors> cpnGraph; //compensation graph
        unsigned begin_idx,lst_idx;
    };

    struct MergeType{
    public:
        std::vector<NNList> knnGraph;
        std::vector<Neighbors> cpnGraph; //compensation graph
    };
}









#endif // PVETEX_H
