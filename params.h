#ifndef PARAMS_H
#define PARAMS_H
#include "string"

struct Params
{
public:
    unsigned S; ///num. of neighbors to be visited
    unsigned init; ///0 by default
    unsigned T; ///num. of iterations
    unsigned P; ///num. of samples
    unsigned K; ///top-K to be kept
};

struct GraphInfo
{
public:


    unsigned N;
    unsigned K;
    unsigned *toplst;
    float    *dstlst;
    unsigned *absId;
    unsigned *relId;
    public:
    GraphInfo():N(0),K(0){}

    GraphInfo(int n, int k):N(n),K(k) {
        dstlst = new float[n*k];
        toplst = new unsigned[n*k];
        absId  = new unsigned[n];
    }

    ~GraphInfo(){

        if(dstlst != NULL)
        {
            delete [] dstlst;
            dstlst = NULL;
        }
        if(toplst != NULL)
        {
            delete [] toplst;
            toplst = NULL;
        }
        if(absId != NULL)
        {
            delete [] absId;
            absId = NULL;
        }
    }
};

#endif
