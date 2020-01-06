#ifndef SPACE_LP_H
#define SPACE_LP_H

#include "abstractgraph.h"
#include "vmath.h"
#include "distopt.h"


 class L2space: public SpaceInterface<float>
    {
        DISTFUNC<float> distfunc;

    public:
        L2space()
        {
            distfunc = VMath::l2norm;
        }

        virtual DISTFUNC<float> get_dist_func()
        {
            return  distfunc;
        }
    };

    class L1space: public SpaceInterface<float>
    {
        DISTFUNC<float> distfunc;

    public:
        L1space()
        {
            distfunc = DistOpt::l1norm;
        }
         virtual DISTFUNC<float> get_dist_func()
        {
            return  distfunc;
        }

    };

    class Ksquare: public SpaceInterface<float>
    {
        DISTFUNC<float> distfunc;
    public:
        Ksquare()
        {
            distfunc = DistOpt::KSquare;
        }
         virtual DISTFUNC<float> get_dist_func()
        {
            return  distfunc;
        }
    };


#endif // SPACE_LP_H
