#ifndef VMATH_H
#define VMATH_H

#include <vector>
#include <random>
#include <set>
#include <map>
#define LG(a,b) (log(a)/log(b))


using namespace std;

class VMath {
private:
    static const float smallVal0;
public:
    static double cos(const float v1[], const unsigned int idx1, const float v2[],
                     const unsigned int idx2, const unsigned int dim0);
    static float  cos(const float v1[], const float v2[], const unsigned int d0);

    static float  l2(const float v1[], const float v2[], const unsigned int d0);

    static double l2(const float v1[], const unsigned int s1, const float v2[],
                    const unsigned int s2, const unsigned int d0);

    static float l2norm(const float *v1, const unsigned long s1,
                        const float *v2, const unsigned long s2, const unsigned int d0);

    static int  sign(const float val) {if(val >= 0) return 1; else return -1; }

    static void randSeed(set<unsigned int> &seeds, unsigned int bound, unsigned int sz);
    static void GenRandom (std::mt19937 &rng, unsigned *addr, unsigned size, unsigned N);

    static int   overlap(const set<int>* seta, const set<int>* setb);
    static unsigned int overlap(const set<unsigned int> &seta, const set<unsigned int> &setb);
    static float normVect(multimap<unsigned int, float> &bowVect);
    static float normVect(float  *feat,  const unsigned int dim);
    static float normVects(float *feats, const unsigned int dim, const unsigned int numb);

    static float lgx(const float a, const float b);
    static float absx(const float val);
    static float inner_product(const float* vect1, const unsigned int size);
    static bool  isZeros(const float *feat, const unsigned int dim,
                         const unsigned int i);
    static bool  isZeros(const float *feat, const unsigned int dim);
    static void  test();
};
#endif
