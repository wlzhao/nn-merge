#ifndef DISTOPT_H
#define DISTOPT_H

#include <unordered_set>

using namespace std;

class DistOpt {
   public:
    static float L2SqrSIMDExt(const float* v1, unsigned int loc1,
                              const float* v2, unsigned int loc2, unsigned int nDim, float tmp[4]);
    static float l1norm(const float *v1, const unsigned long s1, const float *v2,
                          const unsigned long s2, const unsigned int d0);

    static float KSquare(const float* vc1, const unsigned long loc1,
                            const float* vc2, const unsigned long loc2, unsigned int nDim);
    static float cosine(const float *v1, const unsigned int s1, const float *v2,
                          const unsigned int s2, const unsigned int d0);

    static float jaccard(unordered_set<unsigned int> &v1, unordered_set<unsigned int> &v2);

    static void test();
};


#endif // DISTOPT_H_INCLUDED
