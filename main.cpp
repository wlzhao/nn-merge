#include <iostream>

#include "lynndescent.h"
#include "iodelegator.h"
#include "evaluator.h"
#include "lysearch.h"
#include "kgraph.h"
#include "abstractmerge.h"
#include "smerge.h"
#include "jmerge.h"

using namespace std;
using namespace kgraph;

/*****
@author: Peng-Cheng Lin, Wan-Lei Zhao
@date:   Jan.-6-2020

This project is an implementation of "On the Merge of k-NN graph"
that is proposed by Wan-Lei Zhao. Part of the codes are collected
from KGraph project that is maintained by Dr. Wei Dong. Great Thanks
to Dr. Wei Dong!

***/


void test()
{
    Smerge::test(); ///Symmetric Merge
    ///Jmerge::test(); ///Joint Merge
    ///LyNNDescent::test(); ///Hierarchical k -NN graph Construction via J-Merge
    ///LySearch::test();  /// NN search task based on H-Merge
}


int main()
{

    test();
    return 0;
}
