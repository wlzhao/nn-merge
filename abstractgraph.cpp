#include "abstractgraph.h"
#include <algorithm>

float AbstractGraph::default_delta = 0.002;


void AbstractGraph::getRndSds(std::mt19937 &rng, const unsigned int k, unsigned int *addr, const unsigned int N)
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

AbstractGraph::AbstractGraph(SpaceInterface<float> *func)
{
    distfunc = func->get_dist_func();
}


AbstractGraph::~AbstractGraph()
{
    delete this->mat;
    this->mat = NULL;
}
