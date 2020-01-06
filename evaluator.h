#ifndef EVALUATOR_H
#define EVALUATOR_H


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>

using namespace std;

class Evaluator
{
    private:
        static const unsigned int dvd0;
    public:
        Evaluator(){}
        virtual ~Evaluator(){}

        static float NNRecall(const char* srcfn1, const char *grdfn, ostream &out_strm, unsigned int topk0);
        static void NNFailed(const char* srcfn1, const char *grdfn, const char *failedfn, unsigned int topk0);
        static void eval_Clust(const char *clustfn, const char *nmapfn,
                               const char *gtfn, const char *dstfn);

        static float getKGraphRecall(const char *srcFn, const char *grdFn, const unsigned int atop);
        static float getRGraphRecall(const char *srcFn, const char *grdFn);

        static void topk_View(const char *srcfn, const char *grdfn, const char *itms,
                           const char *dstfn, const unsigned int topk0);

        static int  intersect(set<unsigned int> &set1, set<unsigned int> &set2);
        static int  intersect(vector<unsigned int> &vect1, set<unsigned int> &set1, const unsigned int k);
        static int  intersect(vector<unsigned int> &vect1, set<unsigned int> &set1);
        static void loadGrdKNNGraph(const char *srcFn, map<unsigned int, set<unsigned int>* > &kNNGraph, const int topk);
        static void loadGrdRNNGraph(const char *srcFn, map<unsigned int, set<unsigned int>* > &kNNGraph);
        static void loadKNNGraph(const char *srcFn, map<unsigned int, vector<unsigned int>* > &kNNGraph, const int topk);
        static void loadRNNGraph(const char *srcFn, map<unsigned int, vector<unsigned int>* > &kNNGraph);
        static void load_NNmaps(const char* srcfn, map<string, unsigned int> &nmap);
        static void load_NNmap(const char* srcfn, map<unsigned int, unsigned int> &nnmap);
        static void load_NNmap(const char* srcfn, map<string, unsigned int> &nnmap);
        static unsigned int load_clusts(const char *srcfn, map<string, unsigned int> &nmap,
                            map<unsigned int, unsigned int> &grd_clust);
        static void NNRecall(const char *srcfn, const char *groundtruth,
                             const unsigned int topk0, const char *dstfn);
        static void clear_Clusts(map<unsigned int, set<unsigned int> *> &clusts);
        static void clear_grdmap(map<string, set<string>* > &grdmap);
        static void clear_NameMap(map<string, unsigned int> &nmap);
        static void test();
        static void testYJ();
};

#endif
