#include "vmath.h"

#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cstdio>

#include <cstring>
#include <cmath>


using namespace std;

const float VMath::smallVal0 = 0.00005;

float VMath::l2(const float v1[], const float v2[], const unsigned int d0)
{
    assert(v1);
    assert(v2);

    float w1  = 0.0f, w2 = 0.0f;
    float val = 0.0f;

    for(unsigned int i = 0; i < d0; i++)
    {
        w1 = v1[i] - v2[i];
        w2 += w1*w1;
    }
    val = sqrt(w2);
    return val;
}

double VMath::l2(const float v1[], unsigned int s1, const float v2[],
                 unsigned int s2,  const unsigned int d0)
{
    assert(v1);
    assert(v2);

    double w1  = 0.0f, w2 = 0.0f;
    double val = 0.0f;
    unsigned int loc1 = d0 * s1, i = 0;
    unsigned int loc2 = d0 * s2;

    for(i = 0; i < d0; i++)
    {
        w1 = v1[loc1+i] - v2[loc2+i];
        w2 += w1*w1;
    }

    val = w2;

    return val;
}

float VMath::l2norm(const float *v1, const unsigned long s1,
                    const float *v2, const unsigned long s2, const unsigned int d0)
{
    float res = 0, diff = 0;
    const float *pVect1 = v1 + s1*d0;
    const float *pVect2 = v2 + s2*d0;

    size_t qty4 = d0/4;
    const float* pEnd1 = pVect1 + qty4 * 4;
    const float* pEnd2 = pVect1 + d0;

    while (pVect1 < pEnd1)
    {
        diff = *pVect1++ - *pVect2++;
        res += diff * diff;
        diff = *pVect1++ - *pVect2++;
        res += diff * diff;
        diff = *pVect1++ - *pVect2++;
        res += diff * diff;
        diff = *pVect1++ - *pVect2++;
        res += diff * diff;
    }



    while (pVect1 < pEnd2)
    {
        diff = *pVect1++ - *pVect2++;
        res += diff * diff;
    }

    return res;
}

double VMath::cos(const float v1[], const unsigned int idx1, const float v2[],
                  const unsigned int idx2, const unsigned int d0)
{
    assert(v1);
    assert(v2);

    double w1  = 0.0f, w2 = 0.0f;
    double val = 0.0f;
    unsigned int loc1 = d0 * idx1;
    unsigned int loc2 = d0 * idx2;

    for(unsigned int i = 0; i < d0; i++)
    {
        w1 = w1 + v1[loc1+i]*v1[loc1+i];
        w2 = w2 + v2[loc2+i]*v2[loc2+i];
        val = val + v1[loc1+i]*v2[loc2+i];
    }

    w1 = sqrt(w1);
    w2 = sqrt(w2);

    if(w1 == 0 && w2 == 0)
    {
        return 1;
    }

    if(w1 == 0)
    {
        return 0;
    }
    else if(w2 == 0)
    {
        return 0;
    }
    else
    {
        val = val/(w1*w2);
        return val;
    }
}

void VMath::GenRandom (std::mt19937 &rng, unsigned *addr, unsigned size, unsigned N)
{
    if (N == size)
    {
        for (unsigned i = 0; i < size; ++i)
        {
            addr[i] = i;
        }
        return;
    }
    for (unsigned i = 0; i < size; ++i)
    {
        addr[i] = rng() % (N - size);
    }
    sort(addr, addr + size);
    for (unsigned i = 1; i < size; ++i)
    {
        if (addr[i] <= addr[i-1])
        {
            addr[i] = addr[i-1] + 1;
        }
    }
    unsigned off = rng() % N;
    for (unsigned i = 0; i < size; ++i)
    {
        addr[i] = (addr[i] + off) % N;
    }
}


float VMath::cos(const float v1[], const float v2[], const unsigned int d0)
{
    assert(v1);
    assert(v2);

    float w1  = 0.0f, w2 = 0.0f;
    float val = 0.0f;
    for(unsigned int i = 0; i < d0; i++)
    {
        w1 = w1 + v1[i]*v1[i];
        w2 = w2 + v2[i]*v2[i];
        val = val + v1[i]*v2[i];
    }

    w1 = sqrt(w1);
    w2 = sqrt(w2);

    if(w1 == 0 && w2 == 0)
    {
        return 1;
    }

    if(w1 == 0)
    {
        return 0;
    }
    else if(w2 == 0)
    {
        return 0;
    }
    else
    {
        val = val/(w1*w2);
        return val;
    }
}

void VMath::randSeed(set<unsigned int> &seeds, unsigned int bound, unsigned int sz)
{
    unsigned int r, cnt = 0;
    unsigned char *visited = new unsigned char[bound];
    memset(visited, 0, bound);
	do{
	   r = rand()%bound;
	   if(visited[r] == 0)
	   {
	       seeds.insert(r);
	       visited[r] = 1;
	       cnt++;
       }
	}while(cnt < sz);
	delete [] visited;
	visited = NULL;
}


float VMath::inner_product(const float* vect1, const unsigned int size)
{
    assert(vect1);

    float dist = 0;
    for(unsigned int i = 0; i < size; i++)
    {
        dist += vect1[i]*vect1[i];
    }
    return dist;
}

inline float VMath::lgx(const float a,const float b)
{
    return log(a)/log(b);
}

float VMath::absx(const float val)
{
    float  tmpval = val<0?(-1*val):val;
    return tmpval;
}


int VMath::overlap(const set<int> *seta, const set<int> *setb)
{
    assert(seta);
    assert(setb);

    int i = 0;
    set<int>::iterator it;
    for(it = seta->begin(); it != seta->end(); it++)
    {
        if(setb->find(*it) != setb->end())
        {
            i++;
        }
    }
    return i;
}

unsigned int VMath::overlap(const set<unsigned int> &seta,const set<unsigned int> &setb)
{
    unsigned int i = 0;
    set<unsigned int>::iterator it;
    for(it = seta.begin(); it != seta.end(); it++)
    {
        if(setb.find(*it) != setb.end())
        {
            i++;
        }
    }
    return i;
}

float VMath::normVect(multimap<unsigned int, float> &vect)
{
    map<unsigned int, float>::iterator mit;
    float sum = 0;

    for(mit = vect.begin(); mit != vect.end(); mit++)
    {
        sum += mit->second*mit->second;
    }
    sum = sqrtf(sum);
    if(sum > smallVal0)
    {
        for(mit = vect.begin(); mit != vect.end(); mit++)
        {
            mit->second = mit->second/sum;
        }
    }
    return 1.0;
}

float VMath::normVect(float *feats, const unsigned int dim)
{
    assert(feats);
    unsigned int j;
    float len = 0;
    for(j = 0; j < dim; j++)
    {
        len += feats[j]*feats[j];
    }

    len = sqrt(len);
    if(len > 0)
    {
        for(j = 0; j < dim; j++)
        {
            feats[j] = feats[j]/len;
        }
    }
    else
    {
        len = 0.0f;
    }
    return len;
}


float VMath::normVects(float *feats, const unsigned int dim,
                       const unsigned int numb)
{
    assert(feats);
    unsigned int i, j, loc;
    float len = 0;
    for(i = 0; i < numb; i++)
    {
        len = 0;
        loc = i*dim;
        for(j = 0; j < dim; j++)
        {
            len += feats[loc+j]*feats[loc+j];
        }

        len = sqrt(len);
        if(len > smallVal0)
        {
            for(j = 0; j < dim; j++)
            {
                feats[loc+j] = feats[loc+j]/len;
            }
        }
    }
    return 1.0;
}

bool VMath::isZeros(const float *feat, const unsigned int dim,
                    const unsigned int i)
{
    assert(feat);
    unsigned int j, loc;
    float len = 0;
    loc = i*dim;
    for(j = 0; j < dim; j++)
    {
        len += feat[loc+j]*feat[loc+j];
    }

    len = sqrt(len);
    if(len < smallVal0)
    {
        return true;
    }
    else
        return false;
}

bool VMath::isZeros(const float *feat, const  unsigned int dim)
{
    assert(feat);
    unsigned int j;
    float len = 0;
    for(j = 0; j < dim; j++)
    {
        len += feat[j]*feat[j];
    }

    len = sqrt(len);
    if(len < smallVal0)
    {
        return true;
    }
    else
        return false;
}

