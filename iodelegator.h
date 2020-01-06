#ifndef IODELEGATOR_H
#define IODELEGATOR_H

#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <vector>
#include <string>
#include <list>
#include <map>
#include <set>

//#include "invtitem.h"

using namespace std;


class IODelegator
{
private:
    static const unsigned int LONGCHAR  = 2000;
    static const unsigned int FNLEN     = 1024;
    static const unsigned int precision = 3;

public:
    IODelegator();
    static long   getFileSize(const char *srcfn);
    static float *loadMatrix(const char *fn, unsigned int &num, unsigned int &dim);
    static float *loadVocab(const char *fn,  unsigned int &row, unsigned int &col);
    static float *loadXQVocab(const char *vocabFn,  unsigned int &nqDim, unsigned int numb[16], unsigned int &dv0);
    static void   loadMultiMap(const char *srcfn, vector<vector<unsigned int> > &supmap);
    static void   load_NNmap(const char* srcFn, map<unsigned int, unsigned int> &nnmap);

    static void   loadSeeds(const char *srcFn, vector<unsigned int> &mySeeds);

    static float *load_refSet(const char *srcFn,   unsigned int &d, unsigned int &r);
    static float *load_fvecs(const char *fvecFn,   unsigned int &d, unsigned int &r);
    static float *load_bvecs(const char *bvecFn,   unsigned int &d, unsigned int &r);
    static float *loadBvec10m(const char *bvecfn,  unsigned int &d, unsigned int &r);
    static float *loadBvec100m(const char *bvecfn, unsigned int &d, unsigned int &r);
    static void loadSMat(const char *srcFn, unsigned int &r, vector<unordered_set<unsigned> > &smat);

    static float *bvec2txt(const char *bvecfn, const char *dstfn, unsigned int &d, unsigned int &r);

    static unsigned int  *loadKNNGraph(const char *srcFn, unsigned int &gRow, unsigned int &gDim, const unsigned int topk0);
    static unsigned int **loadANNGraph(const char *srcFn, unsigned int &gRow, unsigned int &gDim);

    static void extrTrucMat(const char *srcTrc,      const char *srcFn, const char *dstFn);
    static void extrTrucMat(const unsigned int numb, const char *srcFn, const char *dstFn);
    static float *loadTruncMat(const char *srcFn, unsigned int &row, unsigned int &dim, vector<unsigned int> &ids);
    static void saveTruncMat(const char *dstFn, const unsigned int dim, float *mat, set<unsigned int> &ids);

    static void saveTruncKnnGraph(const char *knnFn, const char *dstFn, unsigned int *idx, unsigned int idxNum);
    static void randSelect(vector<unsigned int> &pool, set<unsigned int> &sel, const unsigned int sz);
    static void procNUSWide(const char *srcFn, const char *qryFn, const char *refFn, const char *trnFn);

    static bool loadWghs(const char *weightfn, map<unsigned int,float> &i2wMap, const unsigned int offset);
    static bool loadWghs(const char *weightfn, map<unsigned int,float> &i2wMap, map<unsigned int, unsigned int>  &pmap,
                         map<unsigned int, unsigned int>  &_pmap, map<unsigned int, int>  &tmap, const unsigned int offset);

    static void load_k2imap(const char *srcfn,  map<string, unsigned int> &referTab);
    static void load_k2imap(const char *srcfn,  map<string, int> &referTab);
    static void load_i2kmap(const char *nmapfn, map<unsigned int, const char*> &i2kMap,
                            const unsigned int offset);
    static void load_word2ref(const char *srcFn, unsigned int w2fMap[128]);

    static void load_map(const char *srcfn, map<string, unsigned char> &referTab);

    static vector<const char*> loadNameList(const char *kflstfn, const char *msg);
    static vector<const char*> loadStrings(const char *kflstfn);

    static unsigned int *loadDPG(const char *dpgFn, unsigned int &row, unsigned int &dim);
    static unsigned int *loadDPGraph(const char *dpgFn, unsigned int &row, unsigned int *&knnLocs);
    static unsigned int *loadTrucDPGraph(const char *dpgFn, unsigned int &row, unsigned int *&knnLocs,
                                         vector<unsigned int> &ids);
    static void test();
    static void testYJ();

    ~IODelegator();
};

#endif
