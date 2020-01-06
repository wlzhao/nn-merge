#ifndef SMERGE_H
#define SMERGE_H
#include "abstractmerge.h"
#include "evaluator.h"



class Smerge: public AbstractMerge
{
    protected:
        float *f_mat, *b_mat;
        unsigned *GraphOne, *GraphTwo;
    public:
        Smerge();
        Smerge(SpaceInterface<float> *func);
        ~Smerge();
        virtual unsigned InitialGraph();
        virtual unsigned* buildKGraph(float *mata, const unsigned K, const unsigned N, const unsigned D);
        virtual unsigned  localJoin();
        virtual void setdisfunc(SpaceInterface<float> *func);

        static void test();
};



#endif // EMERGE_H
