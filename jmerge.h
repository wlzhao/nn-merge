#ifndef JMERGE_H
#define JMERGE_H
#include "abstractmerge.h"

#include "evaluator.h"

class Jmerge: public AbstractMerge
{

    protected:
        float *f_mat;
        unsigned *GraphOne;
    public:
        Jmerge();
        Jmerge(SpaceInterface<float> *func);
        ~Jmerge();
        virtual unsigned InitialGraph();
        virtual unsigned* buildKGraph(float *mata, const unsigned K, const unsigned N, const unsigned D);
        virtual unsigned localJoin();
        virtual void setdisfunc(SpaceInterface<float> *func);
        static void test();
};
#endif // JMERGE_H
