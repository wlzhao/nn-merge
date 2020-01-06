#ifndef CLEANER_H
#define CLEANER_H

#include <vector>
#include <map>
#include <set>

#include "nnitem.h"

using namespace std;

class Cleaner
{
    public:
        Cleaner() {}
        virtual ~Cleaner(){}
        static void clearRNNGraph(map<unsigned int, vector<unsigned int>* > &rnnGraph);
        static void clearRNNGraph(map<unsigned int, vector<NNItem*>* > &rnnGraph);
        static void clearVNNGraph(unsigned int **vnnGraph, const unsigned int nRow);
        static void clearKNNGraph(map<unsigned int, vector<unsigned int>* > &kNNGraph);
        static void clearKNNGraph(map<unsigned int, set<unsigned int>* > &kNNGraph);
};

#endif
