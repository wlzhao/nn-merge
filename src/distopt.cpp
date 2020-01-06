#include "distopt.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <xmmintrin.h>
#include <ammintrin.h>

#include <iostream>
using namespace std;

float DistOpt::L2SqrSIMDExt(const float* vc1, unsigned int loc1,
                            const float* vc2, unsigned int loc2, unsigned int nDim, float TmpRes[4])
{

    size_t qty4 = nDim >> 2;
    size_t qty16 = nDim >> 4;
    const float *pVect1 = vc1 + loc1*nDim;
    const float *pVect2 = vc2 + loc2*nDim;
    const float* pEnd1 = pVect1 + (qty16 << 4);
    const float* pEnd2 = pVect1 + (qty4 << 2);
    const float* pEnd3 = pVect1 + nDim;

    float delta, res;

    __m128  diff, v1, v2;
    __m128  sum = _mm_set1_ps(0);

    while (pVect1 < pEnd1)
    {
        //_mm_prefetch((char*)(pVect2 + 16), _MM_HINT_T0);
        v1 = _mm_loadu_ps(pVect1);
        pVect1 += 4;
        v2 = _mm_loadu_ps(pVect2);
        pVect2 += 4;
        diff = _mm_sub_ps(v1, v2);
        sum = _mm_add_ps(sum, _mm_mul_ps(diff, diff));

        v1 = _mm_loadu_ps(pVect1);
        pVect1 += 4;
        v2 = _mm_loadu_ps(pVect2);
        pVect2 += 4;
        diff = _mm_sub_ps(v1, v2);
        sum = _mm_add_ps(sum, _mm_mul_ps(diff, diff));

        v1 = _mm_loadu_ps(pVect1);
        pVect1 += 4;
        v2 = _mm_loadu_ps(pVect2);
        pVect2 += 4;
        diff = _mm_sub_ps(v1, v2);
        sum = _mm_add_ps(sum, _mm_mul_ps(diff, diff));

        v1 = _mm_loadu_ps(pVect1);
        pVect1 += 4;
        v2 = _mm_loadu_ps(pVect2);
        pVect2 += 4;
        diff = _mm_sub_ps(v1, v2);
        sum = _mm_add_ps(sum, _mm_mul_ps(diff, diff));
    }

    while (pVect1 < pEnd2)
    {
        v1 = _mm_loadu_ps(pVect1);
        pVect1 += 4;
        v2 = _mm_loadu_ps(pVect2);
        pVect2 += 4;
        diff = _mm_sub_ps(v1, v2);
        sum = _mm_add_ps(sum, _mm_mul_ps(diff, diff));
    }

    _mm_store_ps(TmpRes, sum);
    res = TmpRes[0] + TmpRes[1] + TmpRes[2] + TmpRes[3];

    while (pVect1 < pEnd3)
    {
        delta = *pVect1++ - *pVect2++;
        res += delta * delta;
    }

    return res;
};

float DistOpt::l1norm(const float *v1, const unsigned long s1, const float *v2,
                          const unsigned long s2, const unsigned int d0)
{
    float res = 0, diff = 0;
    const float *pVect1 = v1 + s1*d0;
    const float *pVect2 = v2 + s2*d0;

    size_t qty4 = d0/4;
    const float* pEnd1 = pVect1 + qty4 * 4;
    const float* pEnd2 = pVect1 + d0;

    while(pVect1 < pEnd1)
    {
        diff = *pVect1++ - *pVect2++;
        res += fabs(diff);
        diff = *pVect1++ - *pVect2++;
        res += fabs(diff);
        diff = *pVect1++ - *pVect2++;
        res += fabs(diff);
        diff = *pVect1++ - *pVect2++;
        res += fabs(diff);
    }

    while (pVect1 < pEnd2)
    {
        diff = *pVect1++ - *pVect2++;
        res += fabs(diff);
    }

    ///res = sqrt(res);

    return res;
}


float DistOpt::cosine(const float *v1, const unsigned int s1, const float *v2,
                          const unsigned int s2, const unsigned int d0)
{
    float res = 0, l2 = 0, l1 = 0;
    const float *pVect1 = v1 + s1*d0;
    const float *pVect2 = v2 + s2*d0;

    size_t qty4 = d0/4;
    const float* pEnd1 = pVect1 + qty4 * 4;
    const float* pEnd2 = pVect1 + d0;

    while(pVect1 < pEnd1)
    {
        l1  += (*pVect1)*(*pVect1);
        l2  += (*pVect2)*(*pVect2);
        res += (*pVect1++)*(*pVect2++);

        l1  += (*pVect1)*(*pVect1);
        l2  += (*pVect2)*(*pVect2);
        res += (*pVect1++)*(*pVect2++);

        l1  += (*pVect1)*(*pVect1);
        l2  += (*pVect2)*(*pVect2);
        res += (*pVect1++)*(*pVect2++);

        l1  += (*pVect1)*(*pVect1);
        l2  += (*pVect2)*(*pVect2);
        res += (*pVect1++)*(*pVect2++);
    }

    while (pVect1 < pEnd2)
    {
        l1  += (*pVect1)*(*pVect1);
        l2  += (*pVect2)*(*pVect2);
        res += (*pVect1++)*(*pVect2++);
    }

    if(l1 != 0 && l2 != 0)
      res = res/(sqrt(l1)*sqrt(l2));
    else{
       res = 0;
    }
    res = 2 - 2*res;

    return res;
}


float DistOpt::KSquare(const float* v1, const unsigned long s1,
                            const float* v2, const unsigned long s2, unsigned int d0)
{
    float w1  = 0.0f, w2 = 0.0f;
    unsigned int loc1 = d0 * s1, i = 0;
    unsigned int loc2 = d0 * s2;

    for(i = 0; i < d0; i++)
    {
        w1 = v1[loc1+i] - v2[loc2+i];

        if(v1[loc1+i] != 0 || v2[loc2+i] != 0)
        {
            w2 += w1*w1/(v1[loc1+i] + v2[loc2+i]);
        }
        ///cout<<w1<<"\t"<<w2<<"\t"<<v1[loc1+i]<<"\t"<<v2[loc2+i]<<endl;
    }

   return w2;
}


float DistOpt::jaccard(unordered_set<unsigned int> &v1, unordered_set<unsigned int> &v2)
{
    unordered_set<unsigned int>::iterator sit;
    float w = 0, d = 0;
    d = v1.size() + v2.size();
    if(v1.size() == 0 && v2.size() == 0)
    return 0;

    for(sit = v1.begin(); sit != v1.end(); sit++)
    {
        if(v2.find(*sit) != v2.end())
        {
           w = w + 1.0;
        }
    }

    w = w/(d - w);

    return w;
}

void DistOpt::test()
{
     float v1[6] = {1.2, 4.1, 6.1, 2.2, 4.5, 7.1};
     float v2[6] = {2.2, 3.1, 5.1, 3.2, 5.5, 3.3};
     float tmp[4] = {0};
     float l1 = DistOpt::L2SqrSIMDExt(v1, 0, v2, 0, 6, tmp);
     printf("l2 = %f\n", l1);
}
