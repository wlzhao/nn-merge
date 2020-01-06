#include "lysearch.h"

#include "scriptparser.h"
#include "iodelegator.h"
#include "evaluator.h"

#include <iostream>
#include <fstream>
#include <cstring>
#include <cassert>
#include <map>
#include <queue>
#include <memory>

using namespace std;

const unsigned int LySearch::Sd = 40;
const unsigned int LySearch::L0 = 2; //12 for bSearch
const unsigned int LySearch::maxlevel = 5;

LySearch::LySearch(const char *conf)
{
    this->nLayer = 0;
    map<string, const char*> paras = ScriptParser::getConf(conf);

    this->loadRefSet(paras["refer"]);
    this->loadGraphs(paras["dpg"]);

    ScriptParser::clearParaMap(paras);
}

LySearch::LySearch(const char *Rsrc, const char *Gsrc, SpaceInterface<float> *func)
{
    this->nLayer = 0;
    this->loadRefSet(Rsrc);
    this->loadGraphs(Gsrc);
    ///this->loadsingleGraph(Gsrc);
    distfunc = func->get_dist_func();
}

LySearch::LySearch()
{
    this->nCmps = 0;
    this->mat   = NULL;
    this->crntAugKnnGraph = NULL;
    this->crntKnnLocs = NULL;
}

LySearch::LySearch(SpaceInterface<float> *func)
{
    this->nCmps = 0;
    this->mat   = NULL;
    this->crntAugKnnGraph = NULL;
    this->crntKnnLocs = NULL;
    distfunc = func->get_dist_func();

}

int LySearch::loadRefSet(const char *srcFn)
{
    this->mat  = IODelegator::load_refSet(srcFn, this->col, this->row);
    cout<<"Size of reference ................................ "<<this->row<<"x"<<this->col<<endl;
    assert(row > 0);
    return 0;
}

int LySearch::loadsingleGraph(const char *srcFn)
{
    ifstream is(srcFn, ios::in);
    unsigned int i = 0, j = 0, idx = 0, numb, val = 0;
    unsigned long cnts = 0, k = 0;
    assert(is.is_open());
    is>>row;
    this->crntKnnLocs = new unsigned int[row];
    ///cout<<row<<endl;
    this->crntKnnLocs[0] = 0;
    for (i = 0; i < row && !is.eof(); ++i)
    {
        is>>idx;
        is>>numb;
        cnts = cnts + numb + 1;
        for (j = 0; j < numb; ++j)
        {
           is>>val;
        }
        if(i+1 < row)
        {
            this->crntKnnLocs[i+1] = this->crntKnnLocs[i] + numb+1;
        }
        idx = 0; numb = 0;
    }
    is.close();

    ifstream is1(srcFn, ios::in);
    this->crntAugKnnGraph = new unsigned int[cnts];
    is1>>row;

    for (i = k = 0; i < row&& !is.eof(); ++i)
    {
        is1>>idx;
        is1>>numb;
        this->crntAugKnnGraph [k] = numb; k++;
        for(j = 0; j < numb; ++j, k++)
        {
           is1>>this->crntAugKnnGraph[k];
        }
    }
    is1.close();
    return 0;
}

int LySearch::loadGraphs(const char *srcFn)
{
    unsigned int i = 0, j = 0, idx = 0, numb, val = 0, g = 0, rowNum = 0;
    unsigned int gNum = LySearch::maxlevel, pre_sz = 0;
    unsigned long k = 0;
    vector<unsigned long> cnts;

    this->layers = new LyGraph[gNum];
    this->nLayer = gNum;
    for(i = 0; i < gNum; i++)
    {
        this->layers[i].augKnnGraph = NULL;
        this->layers[i].knnLocs = NULL;
    }

    cnts.resize(gNum);

    ifstream *inStrm = new ifstream(srcFn, ios::in);
    if(!inStrm->is_open())
    {
        cout<<"Open file '"<<srcFn<<"' failed!"<<endl;
        return 0;
    }

    layers[0].upSz = 0;
    for(g = 0; g < gNum; g ++)
    {
        cnts[g] = 0;
        (*inStrm)>>rowNum;
        layers[g].upSz    = pre_sz;
        layers[g].szGraph = rowNum;
        layers[g].knnLocs = new unsigned int[rowNum];
        layers[g].knnLocs[0] = 0;
        for (i = 0; i < rowNum && !inStrm->eof(); i++)
        {
            (*inStrm)>>idx;
            (*inStrm)>>numb;
            cnts[g] = cnts[g] + numb + 1;
            for (j = 0; j < numb; j++)
            {
                (*inStrm)>>val;
            }
            if(i+1 < rowNum)
            {
                layers[g].knnLocs[i+1] = layers[g].knnLocs[i] + numb+1;
            }
            idx = 0;
            numb = 0;
        }
        pre_sz = rowNum;
    }

    inStrm->seekg(0, ios::beg);

    for(g = 0; g < gNum; g++)
    {
        layers[g].augKnnGraph = new unsigned int[cnts[g]];
        (*inStrm)>>rowNum;
        printf("Graph-size of layer %d ............................ %d\n", g, layers[g].szGraph);
        for (i = k = 0; i < rowNum && !(*inStrm).eof(); ++i)
        {
            (*inStrm)>>idx;
            (*inStrm)>>numb;
            layers[g].augKnnGraph[k] = numb;
            k++;
            for(j = 0; j < numb; j++, k++)
            {
                (*inStrm)>>layers[g].augKnnGraph[k];
            }
        }
    }
    inStrm->close();
    cnts.clear();
    return 0;
}

void LySearch::releaseGraphs(LyGraph *myLayers, const unsigned int nLayer)
{
    unsigned int i = 0;

    if(myLayers == NULL)
        return ;

    for(i = 0; i < nLayer; i++)
    {
        delete [] myLayers[i].augKnnGraph;
        delete [] myLayers[i].knnLocs;
        myLayers[i].augKnnGraph = NULL;
        myLayers[i].knnLocs = NULL;
    }
    delete [] myLayers;
    myLayers = NULL;
}

void LySearch::getRndSds(std::mt19937 &rng, const unsigned int k, unsigned int *addr, const unsigned int N)
{
    if (N == k)
    {
        for (unsigned i = 0; i < k; ++i)
        {
            addr[i] = i;
        }
        return;
    }
    for (unsigned i = 0; i < k; ++i)
    {
        addr[i] = rng() % (N - k);
    }
    sort(addr, addr + k);
    for (unsigned i = 1; i < k; ++i)
    {
        if (addr[i] <= addr[i-1])
        {
            addr[i] = addr[i-1] + 1;
        }
    }
    unsigned off = rng() % N;
    for (unsigned i = 0; i < k; ++i)
    {
        addr[i] = (addr[i] + off) % N;
    }
}


float LySearch::nndHNSWSrch(unsigned cur_obj, float cur_dst, unsigned maxlayer, vector<RTPair> &results, const unsigned topk)
{
    unsigned int j = 0,nID = 0, cid = 0;
    unsigned int  sz = 0, loc_top = 0;
    float  dst = 0, lwbnd = 0;
    unsigned long loc = 0, ts = 0;

    priority_queue<pair<float, unsigned>, vector<pair<float, unsigned>>, CompareByFirst> topLst;
    priority_queue<pair<float, unsigned>, vector<pair<float, unsigned>>, CompareByFirst> candis;

    this->crntAugKnnGraph = this->layers[maxlayer-1].augKnnGraph;
    this->crntKnnLocs     = this->layers[maxlayer-1].knnLocs;

    pair<float, unsigned> crntNode;

    long unsigned avgVisit = 0;

    candis.emplace(-cur_dst, cur_obj);
    topLst.emplace(cur_dst, cur_obj);
    flags[cur_obj] = 1;
    pts.push_back(cur_obj);

    lwbnd = topLst.top().first;

    while(!candis.empty())
    {
        crntNode = candis.top();

        if((-crntNode.first) > lwbnd)
        {
            break;
        }

        candis.pop();
        nID = crntNode.second;
        loc = crntKnnLocs[nID];
        sz  = crntAugKnnGraph[loc];

        for(j = 1; j <= sz; j++)
        {
            cid = crntAugKnnGraph[loc + j];
            if (!flags[cid])
            {
                flags[cid] = 1;
                pts.push_back(cid);
                dst = distfunc(crntQry, 0, mat, cid, col);
                ++nCmps;

                if (topLst.top().first > dst || topLst.size() < topk)
                {
                    candis.emplace(-dst, cid);
                    topLst.emplace(dst, cid);

                    if (topLst.size() > topk)
                    {
                        topLst.pop();
                    }

                    lwbnd = topLst.top().first;
                }
            }
        }///for(j)
    }///while(!candis)



    results.resize(topLst.size());
    loc_top = topLst.size() - 1;

    while(topLst.size() > 0 && loc_top >= 0)
    {
        pair<float, unsigned> rez = topLst.top();
        results[loc_top].dst = rez.first;
        results[loc_top].id  = rez.second;
        --loc_top;
        topLst.pop();
    }

    for(auto vit = pts.begin(); vit != pts.end(); vit++)
    {
        flags[*vit] = 0;
    }

    pts.clear();

    return 0;
}


int LySearch::initconfig(const char *refer, const char *Graph, const char *query)
{
    this->mat = IODelegator::load_refSet(refer,this->col,this->row);
    this->loadsingleGraph(Graph);
    this->queries = IODelegator::load_refSet(query, this->qDim, this->nQry);
    return 1;
}

void LySearch::saveBasedlayer(const char *dst)
{

    unsigned maxlayer = LySearch::maxlevel;

    this->crntAugKnnGraph = this->layers[maxlayer-1].augKnnGraph;
    this->crntKnnLocs     = this->layers[maxlayer-1].knnLocs;

    ofstream *outStrm = new ofstream(dst, ios::out);

    (*outStrm)<<this->row<<endl;

    for(unsigned i = 0; i < this->row; i ++)
    {
        unsigned loc = this->crntKnnLocs[i];
        unsigned maxM = this->crntAugKnnGraph[loc];
        (*outStrm)<<i<<" "<<maxM<<" ";
        for(unsigned j = 1; j <= maxM; j ++)
        {
            (*outStrm)<<this->crntAugKnnGraph[loc + j]<<" ";
        }
        (*outStrm)<<endl;
    }

    outStrm->close();
    delete outStrm;

}


int  LySearch::knnSearch(const char *qryFn, const char *dstFn,const char *grdFn, const char *record)
{
    unsigned i = 0, topk = 0, loc = 0, maxM = 0;
    unsigned l = 0, ir = 0, nn   = 0, maxlayer = 0, t = 0;
    unsigned cur_obj =0, enterpoint = 0;
    long unsigned subtime = 0;
    float dst = 0.0f, cur_dst = 0, r = 0;
    vector<RTPair> topLst;
    float *queries = IODelegator::load_refSet(qryFn, this->qDim, this->nQry);
    ofstream *out = new ofstream(record, ios::out);
    assert(qDim == this->col);
    maxlayer = LySearch::maxlevel;
    cout<<"Start k-NN Graph nn Search ....... \n";
    flags = new unsigned short int[this->row];
    memset(flags, 0, sizeof(unsigned short int)*this->row);


    unsigned Ks1[25] = {60,50,40,36,32,28,24,20, 18 ,15, 14, 12, 10,8,6,4,2};

    enterpoint = rand() % this->layers[0].szGraph;

    for(t = 0; t < 17; t += 1)
    {
        ofstream *outStrm = new ofstream(dstFn, ios::out);
        topk    = Ks1[t];
        subtime = 0;
        nCmps   = 0;
        clock_t t1 = clock();
        for(i = 0; i < nQry; i++ )
        {
            crntQry = queries + i*this->qDim;
            cur_obj = enterpoint;
            cur_dst = distfunc(crntQry, 0, mat, cur_obj, col);
            ++nCmps;
            for(ir = 0; ir <  4; ir ++)
            {
                this->crntAugKnnGraph = this->layers[ir].augKnnGraph;
                this->crntKnnLocs     = this->layers[ir].knnLocs;
                bool flag = true;
                while(flag)
                {
                    loc  = this->crntKnnLocs[cur_obj];
                    maxM = this->crntAugKnnGraph[loc];
                    flag = false;
                    for(l = 1; l <=maxM; l ++)
                    {
                        nn  = crntAugKnnGraph[loc + l];
                        dst = distfunc(crntQry, 0, mat, nn, col);
                        ++nCmps;
                        if(dst < cur_dst)
                        {
                            cur_dst = dst;
                            cur_obj = nn;
                            flag    = true;
                        }
                    }
                }
            }
            nndHNSWSrch(cur_obj, cur_dst, maxlayer, topLst, topk);
            clock_t t2 = clock();
            printTops(i, topLst, outStrm, 1);
            subtime += clock() - t2;
        }
        total_time = (clock() - t1 - subtime)/(CLOCKS_PER_SEC+0.0f);
        r = Evaluator::NNRecall(dstFn, grdFn, cout, 1);
        cout<<this->nCmps/nQry<<"\t"<<this->total_time<<"\n";
        (*out)<<topk<<"\t"<<r<<"\t"<<this->total_time<<"\t"<<this->nCmps/nQry<<endl;
        outStrm->close();
        if(r==1) break;
    }
    delete [] queries;
    queries = NULL;
    ///outStrm->close();
    return 0;
}

float LySearch::nndFlatSrch(unsigned maxlayer, vector<RTPair> &results, const unsigned topk)
{
    unsigned int i = 0, cur_obj = 0;
    unsigned int j = 0, nID = 0, cid = 0;
    unsigned int  sz = 0, loc_top = 0;
    float  dst = 0, lwbnd = 0, cur_dst = 0;
    unsigned long loc = 0;

    priority_queue<pair<float, unsigned>, vector<pair<float, unsigned>>, CompareByFirst> topLst;
    priority_queue<pair<float, unsigned>, vector<pair<float, unsigned>>, CompareByFirst> candis;


    this->crntAugKnnGraph = this->layers[maxlayer-1].augKnnGraph;
    this->crntKnnLocs     = this->layers[maxlayer-1].knnLocs;

    pair<float, unsigned> crntNode;

    mt19937 rng(time(NULL));
    unsigned int *seeds = new unsigned [topk];



    this->getRndSds(rng, topk, seeds, this->row);

    while(topLst.size() < topk)
    {
       cur_obj  = seeds[i];
       cur_dst = distfunc(crntQry, 0, mat, cur_obj, col);
       candis.emplace(-cur_dst, cur_obj);
       topLst.emplace(cur_dst, cur_obj);
       flags[cur_obj] = 1;
       pts.push_back(cur_obj);
       ++i;
       ++nCmps;
    }

    {
        lwbnd = topLst.top().first;
        while(!candis.empty())
        {
            crntNode = candis.top();

            if((-crntNode.first) > lwbnd)
            {
                break;
            }

            candis.pop();
            nID = crntNode.second;
            loc = crntKnnLocs[nID];
            sz  = crntAugKnnGraph[loc];

            for(j = 1; j <= sz; j++)
            {
                cid = crntAugKnnGraph[loc + j];
                if (!flags[cid])
                {
                    flags[cid] = 1;
                    pts.push_back(cid);
                    dst = distfunc(crntQry, 0, mat, cid, col);
                    ++nCmps;

                    if (topLst.top().first > dst || topLst.size() < topk)
                    {
                        candis.emplace(-dst, cid);
                        topLst.emplace(dst, cid);

                        if (topLst.size() > topk)
                        {
                            topLst.pop();
                        }

                        lwbnd = topLst.top().first;
                    }
                }
            }///for(j)
        }///while(!candis)

        results.resize(topLst.size());
        loc_top = topLst.size() - 1;

        while(topLst.size() > 0 && loc_top >= 0)
        {
            pair<float, unsigned> rez = topLst.top();
            results[loc_top].dst = rez.first;
            results[loc_top].id  = rez.second;
            --loc_top;
            topLst.pop();
        }
    }///for(qi)


    for(auto vit = pts.begin(); vit != pts.end(); vit++)
    {
        flags[*vit] = 0;
    }

    pts.clear();

    delete [] seeds;
    seeds = NULL;



    return 0;
}


 int LySearch::singleknnSearch(const char *qryFn, const char *dstFn, const char *grdFn,const char *record)
 {
    unsigned i = 0,topk =0, t = 0, maxlayer = 0;
    long unsigned subtime = 0;
    float r = 0;
    vector<RTPair> topLst;

    float *queries = IODelegator::load_refSet(qryFn, this->qDim, this->nQry);

    ofstream *out = new ofstream(record, ios::out);

    cout<<qDim<<"\t"<<this->col<<endl;
    assert(qDim == this->col);

    maxlayer = LySearch::maxlevel;

    cout<<"Start k-NN Graph nn Search ....... \n";

    unsigned Ks1[25] = {60,50,40,36,32,28,24,20, 18 ,15, 14, 12, 10,8,6,4,2};

    flags = new unsigned short int[this->row];
    memset(flags, 0, sizeof(unsigned short int)*this->row);

    for(t = 0; t < 17; t += 1)
    {
        ofstream *outStrm = new ofstream(dstFn, ios::out);
        topk = Ks1[t];
        subtime = 0;
        nCmps = 0;
        clock_t t1 = clock();
        for(i = 0; i < nQry; i++ )
        {
            crntQry = queries + i*qDim;
            nndFlatSrch(maxlayer, topLst, topk);

            clock_t t2 = clock();
            printTops(i, topLst, outStrm, 1);
            subtime += clock() - t2;
            topLst.clear();
        }

        total_time = (clock() - t1- subtime)/(CLOCKS_PER_SEC+0.0f);
        r = Evaluator::NNRecall(dstFn, grdFn, cout, 1);
        cout<<nCmps/nQry<<"\t"<<total_time<<"\n";
        (*out)<<topk<<"\t"<<r<<"\t"<<total_time<<"\t"<<nCmps/nQry<<endl;
        outStrm->close();
        if(r==1) break;
    }
    delete [] queries;
    queries = NULL;
    return 0;

 }



unsigned LySearch::insert2List(RTPair *addr, unsigned K, RTPair &nwnn)
{
    /// find the location to insert
    int j = 0;
    int i = 0;

    if(nwnn.dst > addr[K - 1].dst) return K + 1;  ///can not be insert.

    j = K;

    while (j > 0)
    {
        j --;
        if (addr[j].dst <= nwnn.dst) break;
        i = j;
    }
    // check for equal ID
    unsigned l = K - 1;
    while (l > 0)
    {
        if (addr[l].dst < nwnn.dst) break;
        if (addr[l].id == nwnn.id) return K + 1;
        l --;
    }
    // i <= K-1
    l = K - 1;

    while (l > i)
    {
        addr[l] = addr[l-1];
        --l;
    }

    addr[i] = nwnn;
    return i;
}

unsigned LySearch::insert2List(RTPair *addr, unsigned K, unsigned id, float dst, const unsigned Sz)
{
    // find the location to insert
    int j = 0;
    unsigned i = 0;

    if(dst > addr[K - 1].dst && K == Sz)
    {
        return K;  ///can not be insert.
    }

    j = K;

    while (j > 0)
    {
        j--;
        if (addr[j].dst < dst) break;
        i = j;
    }

    if(j == K -1)
    {
        i = K;
    }

    int l = 0;

    l = (K == Sz)?(K - 1):K;
    //l = K - 1;

    while (l > i)
    {
        addr[l] = addr[l-1];
        --l;
    }

    addr[i].dst  = dst;
    addr[i].id   = id;
    addr[i].flag = 0;

    return i;
}


void LySearch::printTops(unsigned int i, vector<RTPair> &pairs, ostream *outStrm, int topk)
{
    int k = 0;
    for(k = 0; k < topk; k++)
    {
        (*outStrm)<<i<<"\t"<<pairs[k].id<<"\t"<<pairs[k].dst<<endl;
    }
}

LySearch::~LySearch()
{

    if(this->mat != NULL)
    {
        delete [] this->mat;
        this->mat = NULL;
    }

    releaseGraphs(this->layers, nLayer);
}

void LySearch::savebase()
{
    const char *qryFn  = "../../data/sift_query.txt";
    const char *dstFn  = "../../data/sift_queryResult.txt";
    const char *Rsrc   = "../../data/sift100k.txt";
    const char *Gsrc   = "../../data/sift100k_HMGD_40.txt";
    const char *rndgrd = "../../data/sift_100k_grd.txt";
    const char *record = "../../data/sift_100k_record.txt";


    L2space l2;
    const char *bGsrc   = "../../data/sift_100k_HMGD_baselayer.txt";
    LySearch *mysearch = new LySearch(Rsrc, Gsrc, &l2);
    mysearch->saveBasedlayer(bGsrc);

}


void LySearch::test()
{
    const char *qryFn  = "../../data/sift_query.txt";
    const char *dstFn  = "../../data/sift_queryResult.txt";
    const char *Rsrc   = "../../data/sift100k.txt";
    const char *Gsrc   = "../../data/sift100k_HMGD_40.txt";
    const char *rndgrd = "../../data/sift_100k_grd.txt";
    const char *record = "../../data/sift_100k_record.txt";

    L2space l2;
    ///Ksquare ki;
    LySearch *mysearch = new LySearch(Rsrc, Gsrc, &l2);
    mysearch->knnSearch(qryFn, dstFn, rndgrd,record);
}


void LySearch::flatsearch()
{
    const char *qryFn  = "../../data/sift_query.txt";
    const char *dstFn  = "../../data/sift_queryResult.txt";
    const char *Rsrc   = "../../data/sift100k.txt";
    const char *Gsrc   = "../../data/sift100k_HMGD_40.txt";
    const char *rndgrd = "../../data/sift_100k_grd.txt";
    const char *record = "../../data/sift_100k_record.txt";

    float r = 0;
    L2space l2;
    ///LySearch *mysearch = new LySearch(conf);
    ///Ksquare ki;

    LySearch *mysearch = new LySearch(Rsrc, Gsrc, &l2);
    mysearch->singleknnSearch(qryFn, dstFn, rndgrd,record);

    cout<<"Recall@1: " << r <<endl;

}
