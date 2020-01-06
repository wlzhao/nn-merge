#include "cleaner.h"


using namespace std;

void Cleaner::clearRNNGraph(map<unsigned int, vector<unsigned int>* > &rnnGraph)
{
   vector<unsigned int> *crntVect = NULL;

   for(auto it = rnnGraph.begin(); it != rnnGraph.end(); it++)
   {
         crntVect = it->second;
         crntVect->clear();
         delete crntVect;
   }
   rnnGraph.clear();
}

void Cleaner::clearRNNGraph(map<unsigned int, vector<NNItem*>* > &rnnGraph)
{
   vector<NNItem*>::iterator vit;
   vector<NNItem*> *crntVect = NULL;
   NNItem *crntItm = NULL;

   for(auto it = rnnGraph.begin(); it != rnnGraph.end(); it++)
   {
         crntVect = it->second;
         for(vit = crntVect->begin(); vit != crntVect->end(); vit++)
         {
             crntItm = *vit;
             delete crntItm;
         }
         crntVect->clear();
         delete crntVect;
   }
   rnnGraph.clear();
}

void Cleaner::clearKNNGraph(map<unsigned int, vector<unsigned int>* > &kNNGraph)
{
    map<unsigned int, vector<unsigned int>* >::iterator mit;
    vector<unsigned int>* crntItms = NULL;

    for(mit = kNNGraph.begin(); mit != kNNGraph.end(); mit++)
    {
        crntItms = mit->second;
        if(crntItms != NULL)
        {
          crntItms->clear();
          delete crntItms;
        }
    }
    kNNGraph.clear();
}


void Cleaner::clearKNNGraph(map<unsigned int, set<unsigned int>* > &kNNGraph)
{
    map<unsigned int, set<unsigned int>* >::iterator mit;
    set<unsigned int>* crntItms = NULL;

    for(mit = kNNGraph.begin(); mit != kNNGraph.end(); mit++)
    {
        crntItms = mit->second;
        if(crntItms != NULL)
        {
           crntItms->clear();
           delete crntItms;
        }
    }
    kNNGraph.clear();
}

void Cleaner::clearVNNGraph(unsigned int **vnnGraph, const unsigned int nRow)
{
   unsigned int i = 0;
   unsigned int *nnLst = NULL;
   for(i = 0; i < nRow; i++)
   {
        nnLst = vnnGraph[i];
        if(nnLst != NULL)
        {
           delete [] nnLst;
           nnLst = NULL;
        }
   }
}
