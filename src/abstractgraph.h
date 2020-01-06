#ifndef ABSGRAPH_H
#define ABSGRAPH_H
#include <random>
#include <vector>
#include "vmath.h"
#include "distopt.h"

using namespace std;


template<typename MTYPE>
    using DISTFUNC = MTYPE(*)(const float *, const unsigned long s1, const float *, const unsigned long s2, const unsigned int d0);

    template<typename MTYPE>
    class SpaceInterface {
    public:
        virtual DISTFUNC<MTYPE> get_dist_func() = 0;
    };


class AbstractGraph
{

    protected:
            static float default_delta;
            unsigned K0;
            unsigned row, col;
            float *mat;
            DISTFUNC<float> distfunc;


    public:
            unsigned long nCmps;

            AbstractGraph(): K0(0),row(0), col(0), nCmps(0), mat(NULL) {};
            AbstractGraph(SpaceInterface<float> *func);
            ~AbstractGraph();
            void     getRndSds(std::mt19937 &rng, const unsigned int k, unsigned int *addr, const unsigned int N);
            virtual unsigned InitialMemory(float *mata, const unsigned k, const unsigned n,const unsigned d) = 0;
            virtual unsigned InitialGraph() = 0;
            virtual unsigned* buildKGraph(float *mata, const unsigned K, const unsigned N, const unsigned D) = 0;
            virtual void saveKGraph(const char *dst) = 0;
            virtual void setdisfunc(SpaceInterface<float> *func) = 0;

};







#endif // ABSGRAPH_H
