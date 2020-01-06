#include "evaluator.h"

#include "iodelegator.h"
#include "cleaner.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <cmath>
#include <map>

/**
*
*
1. Evaluate approximate nearest neighbor search approach

2. Evaluate image retrieval performance with mAP
*
*
**/

using namespace std;

const unsigned int Evaluator::dvd0 = 10000;

void Evaluator::NNRecall(const char *srcfn, const char *groundtruth,
                         const unsigned int topk0, const char *dstfn)
{
    unsigned int iline = 0, k;
    unsigned int inn, nn, nquery, nprop, ntopk;
    unsigned int iquery[3];
    float dist;
    float *topk_eval = new float[topk0];

    /**load ground-truth **/
    ifstream *inStrm_grd = new ifstream(groundtruth, ios::in);
    if(!inStrm_grd->is_open())
    {
        cout<<"Grountruth file '"<<groundtruth<<"' cannot be open!\n";
        exit(0);
    }

    (*inStrm_grd)>>nquery;
    (*inStrm_grd)>>nprop;
    (*inStrm_grd)>>ntopk;
    cout<<nquery<<" "<<ntopk<<endl;
    cout<<"Loading groundtruth ...... ";

    map<unsigned int, unsigned char> *topk_map = new map<unsigned int, unsigned char>[nquery];
    map<unsigned int, unsigned char> *crnt_map;

    while(!inStrm_grd->eof() && iline < nquery)
    {
        for(inn = 0; inn < ntopk; inn++)
        {
            (*inStrm_grd)>>nn;
            if(inn < 1)
            {
                topk_map[iline].insert(pair<unsigned int, unsigned char>(nn, 0));
            }/**/
        }
        iline++;
    }
    inStrm_grd->close();
    cout<<iline<<"\n";

    /**parse source file **/
    ifstream *inStrm_src = new ifstream(srcfn, ios::in);
    if(!inStrm_src->is_open())
    {
        cout<<"Source file '"<<srcfn<<"' cannot be open!\n";
        exit(0);
    }

    iline = 0;
    memset(topk_eval, 0, sizeof(float)*topk0);
    while(!inStrm_src->eof())
    {
        crnt_map = &topk_map[iline];
        for(inn = 0; inn < topk0; inn++)
        {
            (*inStrm_src)>>iquery[0];
            (*inStrm_src)>>iquery[1];
            (*inStrm_src)>>iquery[2];
            (*inStrm_src)>>dist;
            if(crnt_map->find(iquery[2]) != crnt_map->end())
            {
                for(k = inn; k < topk0; k++)
                {
                    topk_eval[k] = topk_eval[k] + 1;
                }
            }
        }
        iline++;
    }

    float div = topk0*nquery;
    ofstream *out_strm = new ofstream(dstfn, ios::app);
    if(!out_strm->is_open())
    {
        cout<<"Destine file '"<<dstfn<<"' cannot open!\n";
        exit(0);
    }
    for(inn = 0; inn < topk0; inn++)
    {
        topk_eval[inn] = topk_eval[inn]/div;
        (*out_strm)<<inn+1<<"\t"<<topk_eval[inn]<<endl;
    }
    out_strm->close();

    delete [] topk_eval;
    for(inn = 0; inn < nquery; inn++)
    {
        topk_map[inn].clear();
    }
    delete [] topk_map;
}

void Evaluator::NNFailed(const char* srcfn1, const char *grdfn, const char *failedfn, unsigned int topk0)
{
    map<unsigned int, unsigned int>::iterator mit;
    map<unsigned int, unsigned int> nnmap;
    cout<<"\nLoading ground-truth .................... ";
    Evaluator::load_NNmap(grdfn, nnmap);
    cout<<nnmap.size()<<endl;
    unsigned int qItm[2];
    float dist = 0;
    float recall = 0.0f, prec = 0.0f;
    float total  = (float)nnmap.size();
    ifstream *inStrm = new ifstream(srcfn1, ios::in);
    if(!inStrm->is_open())
    {
       cout<<"File '"<<srcfn1<<"' cannot open for reading!\n";
       exit(0);
    }
    unsigned int counter = 0, nlines = 0;
    map<unsigned int, int> queries;
    bool flag[10000];
    memset(flag, 0, 10000*sizeof(bool));

    ofstream failStrm(failedfn, ios::out);
    if(!failStrm.is_open())
    {
        cout <<"file "<<failedfn << "can't open."<< endl;
        exit(0);
    }

    while(!inStrm->eof())
    {
        (*inStrm)>>qItm[0];
        (*inStrm)>>qItm[1];
        (*inStrm)>>dist;

        if(queries.find(qItm[0]) == queries.end())
        {
            queries.insert(pair<unsigned int, unsigned int>(qItm[0], 1));
            counter = 1;
        }else{
            counter++;
            queries[qItm[0]]++;
        }

        if(queries[qItm[0]] > topk0)
        continue;

        nlines++;

        if(nnmap.find(qItm[0]) != nnmap.end())
        {
            if(qItm[1] == nnmap[qItm[0]])
            {
                recall = recall + 1.0f;
                  prec = prec + 1.0f;
                  flag[qItm[0]] = 1;
            }
        }

    }
    queries.clear();

    inStrm->close();
    cout<<"Top-"<<topk0<<" : "<<recall/total<<"\t"<<nlines/topk0<<endl;

    for(int i = 0; i < 10000; i++)
    {
        if(!flag[i])
        {
            failStrm << i <<" ";
        }
    }
    failStrm.close();

    return;
}

float Evaluator::NNRecall(const char* srcfn1, const char *grdfn, ostream &outStrm, unsigned int topk0)
{
    map<unsigned int, unsigned int>::iterator mit;
    map<unsigned int, unsigned int> nnmap;
    cout<<"\nLoading ground-truth .................... ";
    Evaluator::load_NNmap(grdfn, nnmap);
    cout<<nnmap.size()<<endl;
    unsigned int qItm[2];
    float dist = 0;
    float recall = 0.0f, prec = 0.0f;
    float total  = (float)nnmap.size();
    ifstream *inStrm = new ifstream(srcfn1, ios::in);
    if(!inStrm->is_open())
    {
       cout<<"File '"<<srcfn1<<"' cannot open for reading!\n";
       exit(0);
    }
    unsigned int counter = 0, nlines = 0;
    map<unsigned int, int> queries;

    while(!inStrm->eof())
    {
        (*inStrm)>>qItm[0];
        (*inStrm)>>qItm[1];
        (*inStrm)>>dist;

        if(queries.find(qItm[0]) == queries.end())
        {
            queries.insert(pair<unsigned int, unsigned int>(qItm[0], 1));

            counter = 1;
        }else{
            counter++;
            queries[qItm[0]]++;
        }

        if(queries[qItm[0]] > topk0)
        continue;

        nlines++;

        if(nnmap.find(qItm[0]) != nnmap.end())
        {
            if(qItm[1] == nnmap[qItm[0]])
            {
                recall = recall + 1.0f;
                  prec = prec + 1.0f;
                ///cout<<cnt++<<"\t"<<qItm[0]<<"\t"<<qItm[1]<<endl;
            }
        }
    }
    queries.clear();

    inStrm->close();
    outStrm<<"Top-"<<topk0<<" : "<<recall/total<<"\t"<<nlines/topk0<<endl;

    return recall/total;
}

float Evaluator::getKGraphRecall(const char *srcFn, const char *grdFn, const unsigned int atop)
{
    float r = 0;
    unsigned int key = 0, counts = 0, nhits = 0, hit = 0;
    map<unsigned int, vector<unsigned int>* > kNNGraph;
    map<unsigned int, set<unsigned int>* > grdNNGraph;
    map<unsigned int, vector<unsigned int>* >::iterator mit;
    vector<unsigned int>* crntVect;
    set<unsigned int>* crntSet;
    set<unsigned int>::iterator sit;
    vector<unsigned int>::iterator vit;

    Evaluator::loadKNNGraph(srcFn, kNNGraph, atop);
    Evaluator::loadGrdKNNGraph(grdFn, grdNNGraph, atop);

    for(mit = kNNGraph.begin(); mit != kNNGraph.end(); mit++)
    {
        key = mit->first;
        crntVect = mit->second;
        crntSet  = grdNNGraph[key];

        if(crntSet != NULL)
        {
             hit = Evaluator::intersect(*crntVect, *crntSet, atop);
             nhits += hit;
             counts++;
        }
    }
    r  = nhits/(counts*atop+0.0);

    cout<<"Recall@"<<atop<<" is: "<<r<<endl;

    Cleaner::clearKNNGraph(grdNNGraph);
    Cleaner::clearKNNGraph(kNNGraph);


    return r;
}

float Evaluator::getRGraphRecall(const char *srcFn, const char *grdFn)
{
    float r = 0;
    unsigned int key = 0, hit = 0;
    unsigned long nhits = 0;
    double total = 0;

    map<unsigned int, vector<unsigned int>* > rNNGraph;
    map<unsigned int, set<unsigned int>* > grdNNGraph;
    map<unsigned int, vector<unsigned int>* >::iterator mit;
    vector<unsigned int>* crntVect;
    set<unsigned int>* crntSet;

    Evaluator::loadRNNGraph(srcFn, rNNGraph);
    Evaluator::loadGrdRNNGraph(grdFn, grdNNGraph);

    for(mit = rNNGraph.begin(); mit != rNNGraph.end(); mit++)
    {
        key = mit->first;
        crntVect = mit->second;
        //cout<<key<<"\t"<<crntVect->size()<<endl;

        crntSet  = grdNNGraph[key];

        if(crntSet == NULL)
        continue;

        total   += crntSet->size();

        if(crntSet != NULL && crntVect != NULL)
        {
             hit = Evaluator::intersect(*crntVect, *crntSet);
             nhits += hit;
        }
    }
    r  = nhits/total;
    cout<<"Recall is: "<<r<<"\tnHits: "<<nhits<<endl;

    Cleaner::clearKNNGraph(grdNNGraph);
    Cleaner::clearKNNGraph(rNNGraph);


    return r;
}


void Evaluator::loadGrdKNNGraph(const char *srcFn, map<unsigned int, set<unsigned int>* > &kNNGraph, const int topk)
{
    unsigned int i = 0, id = 0, numb = 0, cid = 0;
    unsigned int sz, it, dim = 0;
    set<unsigned int>* crntItms = NULL;

    kNNGraph.clear();
    ifstream *inStrm = new ifstream(srcFn, ios::in);
    if(!inStrm->is_open())
    {
        cout<<"File '"<<srcFn<<"' cannot open for read!\n";
        exit(0);
    }

    (*inStrm)>>sz;
    (*inStrm)>>dim;
     ///(*inStrm)>>numb;
    for(it = 0; it < sz; it++)
    {
        (*inStrm)>>id;
        (*inStrm)>>numb;
        crntItms = new set<unsigned int>;
        kNNGraph.insert(pair<unsigned int, set<unsigned int>*>(id, crntItms));
        for(i = 0; i < numb; i++)
        {
            (*inStrm)>>cid;
            if(i < topk)
            crntItms->insert(cid);
        }
        ///cout<<id<<"\t"<<crntItms->size()<<endl;
    }

    inStrm->close();

    return ;
}


void Evaluator::loadGrdRNNGraph(const char *srcFn, map<unsigned int, set<unsigned int>* > &kNNGraph)
{
    unsigned int i = 0, id = 0, numb = 0, cid = 0;
    unsigned int sz, it = 0;
    set<unsigned int>* crntItms = NULL;

    kNNGraph.clear();
    ifstream *inStrm = new ifstream(srcFn, ios::in);
    if(!inStrm->is_open())
    {
        cout<<"File '"<<srcFn<<"' cannot open for read!\n";
        exit(0);
    }

    (*inStrm)>>sz;
    for(it = 0; it < sz; it++)
    {
        (*inStrm)>>id;
        (*inStrm)>>numb;
        crntItms = new set<unsigned int>;
        kNNGraph.insert(pair<unsigned int, set<unsigned int>*>(id, crntItms));
        for(i = 0; i < numb; i++)
        {
            (*inStrm)>>cid;
            crntItms->insert(cid);
        }
    }

    inStrm->close();
    cout<<"done\n";

    return ;
}

void Evaluator::loadKNNGraph(const char *srcFn, map<unsigned int, vector<unsigned int>* > &kNNGraph, const int topk)
{
    unsigned int i = 0, id = 0, numb = 0, cid = 0;
    unsigned int sz, it;
    vector<unsigned int>* crntItms = NULL;

    kNNGraph.clear();
    ifstream *inStrm = new ifstream(srcFn, ios::in);
    if(!inStrm->is_open())
    {
        cout<<"File '"<<srcFn<<"' cannot open for read!\n";
        exit(0);
    }

    (*inStrm)>>sz;
    (*inStrm)>>numb;
    ///cout<<"here: "<<sz<<endl;
    for(it = 0; it < sz; it++)
    {
        (*inStrm)>>id;
        (*inStrm)>>numb;

        crntItms = new vector<unsigned int>;
        kNNGraph.insert(pair<unsigned int, vector<unsigned int>*>(id, crntItms));
        for(i = 0; i < numb; i++)
        {
            (*inStrm)>>cid;
            if(i < topk)
            crntItms->push_back(cid);
        }
    }

    inStrm->close();

    return ;
}

void Evaluator::loadRNNGraph(const char *srcFn, map<unsigned int, vector<unsigned int>* > &kNNGraph)
{
    unsigned int i = 0, id = 0, numb = 0, cid = 0;
    unsigned int sz, it;
    vector<unsigned int>* crntItms = NULL;

    kNNGraph.clear();
    ifstream *inStrm = new ifstream(srcFn, ios::in);
    if(!inStrm->is_open())
    {
        cout<<"File '"<<srcFn<<"' cannot open for read!\n";
        exit(0);
    }

    (*inStrm)>>sz;
    (*inStrm)>>numb;
    for(it = 0; it < sz; it++)
    {
        (*inStrm)>>id;
        (*inStrm)>>numb;

        crntItms = new vector<unsigned int>;
        kNNGraph.insert(pair<unsigned int, vector<unsigned int>*>(id, crntItms));
        for(i = 0; i < numb; i++)
        {
            (*inStrm)>>cid;
            crntItms->push_back(cid);
        }
    }

    inStrm->close();
    cout<<"done\n";

    return ;
}

void Evaluator::load_NNmap(const char* srcfn, map<unsigned int, unsigned int> &nnmap)
{
    ifstream inStrm;
    inStrm.open(srcfn, ios::in);
    unsigned int key = 0, val = 0;
    if(!inStrm.is_open())
    {
        cout<<"File '"<<srcfn<<"' cannot open for write!\n";
        return ;
    }

    while(!inStrm.eof())
    {
        inStrm>>key;
        inStrm>>val;
        nnmap.insert(pair<unsigned int, unsigned int>(key, val));
    }
    inStrm.close();
    return ;
}

void Evaluator::load_NNmaps(const char* srcfn, map<string, unsigned int> &nmap)
{
    return ;
}

void Evaluator::load_NNmap(const char* srcfn, map<string, unsigned int> &nnmap)
{
    ifstream inStrm;
    inStrm.open(srcfn, ios::in);
    unsigned int key = 0;
    if(!inStrm.is_open())
    {
        cout<<"File '"<<srcfn<<"' cannot open for write!\n";
        return ;
    }

    char str[512];
    unsigned int offset = nnmap.size();

    while(!inStrm.eof())
    {
        inStrm>>key;
        inStrm>>str;
        string id = str;
        key       = key + offset;
        nnmap.insert(pair<string, unsigned int>(id, key));
    }
    inStrm.close();
    return ;
}

void Evaluator::eval_Clust(const char *clustfn, const char *nmapfn,
                           const char *gtfn, const char *dstfn)
{
    map<unsigned int, unsigned int>  grd_clusts;
    map<string, unsigned int> nmap;

    Evaluator::load_NNmap(nmapfn, nmap);
    Evaluator::load_clusts(gtfn, nmap, grd_clusts);

    grd_clusts.clear();
    Evaluator::clear_NameMap(nmap);
    return ;
}

int Evaluator::intersect(set<unsigned int> &set1, set<unsigned int> &set2)
{
    set<unsigned int>::iterator sit;
    int cmmn = 0;
    for(sit = set1.begin(); sit != set1.end(); sit++)
    {
        if(set2.find(*sit) != set2.end())
        {
            cmmn++;
        }
    }
    return cmmn;
}

int Evaluator::intersect(vector<unsigned int> &vect1, set<unsigned int> &set1, const unsigned int k)
{
     int ovrlap = 0;
     for(auto i = 0; i < k; i++)
     {
         if(set1.find(vect1[i]) != set1.end())
         {
              ovrlap++;
         }
     }
     return ovrlap;
}


int Evaluator::intersect(vector<unsigned int> &vect1, set<unsigned int> &set1)
{
     int ovrlap = 0;
     vector<unsigned int>::iterator vit;
     for(vit = vect1.begin(); vit != vect1.end(); vit++)
     {
         if(set1.find(*vit) != set1.end())
         {
              ovrlap++;
         }
     }
     return ovrlap;
}

unsigned int Evaluator::load_clusts(const char *srcfn, map<string, unsigned int> &nmap,
                            map<unsigned int, unsigned int> &grd_clust)
{
    ifstream *inStrm = new ifstream(srcfn, ios::in);

    cout<<"Loading clusters ........................ ";

    char  str[512], numstr[64];
    int num, i = 0;
    grd_clust.clear();

    assert(inStrm != NULL);

    if(!inStrm->is_open())
    {
        cout<<"File '"<<srcfn<<"' cannot open!\n";
        return 0;
    }
    unsigned int id, clustid = 1;
    unsigned int total_pairs = 0;

    while(!inStrm->eof())
    {
        ///inStrm->getline(line, 8192);
        ///cout<<line<<endl;
        (*inStrm)>>numstr;

        num = atoi(numstr);

        if(num == -1)
        continue;

        total_pairs += num*(num-1)/2;

        for(i = 0; i < num; i++)
        {
            (*inStrm)>>str;

            if(!strcmp(str, ""))
            continue;

            if(nmap.find(str) != nmap.end())
            {
                id = nmap[str];
                if(grd_clust.find(id) == grd_clust.end())
                {
                    grd_clust.insert(pair<unsigned int, unsigned int>(id, clustid));
                }else{
                    cout<<"'"<<str<<"\t"<<clustid<<"' duplicate item!!!!\n";
                }
            }else
            {
                cout<<"Name '"<<str<<"' cannot be found in the map!!!!\n";
            }
        }
        clustid++;
        strcpy(str, "");
        strcpy(numstr, "");
        num = -1;
    }
    inStrm->close();
    cout<<grd_clust.size()<<endl;
    return total_pairs;
}

void Evaluator::clear_grdmap(map<string, set<string>* > &grdmap)
{
    map<string, set<string> *>::iterator mit;
    set<string>::iterator sit;
    set<string> *crnt_lst;

    for(mit = grdmap.begin(); mit != grdmap.end(); mit++)
    {
        string key = mit->first;
        crnt_lst   = mit->second;
        for(sit = crnt_lst->begin(); sit != crnt_lst->end(); sit++)
        {
            string mkey = *sit;
            mkey.clear();
        }
        key.clear();
        crnt_lst->clear();
        delete crnt_lst;
    }
    grdmap.clear();
    return ;
}

void Evaluator::clear_Clusts(map<unsigned int, set<unsigned int> *> &clusts)
{
    map<unsigned int, set<unsigned int> *>::iterator mit;
    set<unsigned int> *crnt_clust;

    for(mit = clusts.begin(); mit != clusts.end(); mit++)
    {
        crnt_clust = mit->second;
        crnt_clust->clear();
        delete crnt_clust;
    }
    clusts.clear();
    return ;
}

void Evaluator::clear_NameMap(map<string, unsigned int> &nmap)
{
    map<string, unsigned int>::iterator mit;
    for(mit = nmap.begin(); mit != nmap.end(); mit++)
    {
        string key = mit->first;
        key.clear();
    }
    nmap.clear();
    return ;
}

void Evaluator::test()
{
    const char *siftgrd  = "/home/wlzhao/datasets/bignn/sift1m/sift_grd128.txt";
    const char *gistgrd  = "/home/wlzhao/datasets/bignn/gist/gist_grd.txt";
    const char *nuswgrd  = "/home/wlzhao/datasets/bignn/nuswide/nusw_grd.txt";
    const char *randgrd  = "/home/wlzhao/datasets/bignn/rand/rand1m100d_grd.txt";
    const char *glovegrd  = "/home/wlzhao/datasets/bignn/glove/glove1m_grd.txt";
    const char *sift100kgrd  = "/home/wlzhao/datasets/bignn/sift1m/sift100k_grd.txt";
    const char *yfccgrd  = "/home/wlzhao/datasets/bignn/yfcc1m/yfcc1m_grd.txt";;
    const char *knngrd1  = "/home/wlzhao/datasets/bignn/sift1m/sift1m_gold_knn30.txt";
    const char *knnsrc1  = "/home/wlzhao/datasets/bignn/sift1m/sift1m_k=40v9.txt";
    const char *knngrd2  = "/home/wlzhao/datasets/bignn/gist1m/gist1m_gold_knn30.txt";
    const char *knnsrc2  = "/home/wlzhao/datasets/bignn/gist1m/gist1m_k=40v9.txt";
    const char *knngrd3  = "/home/wlzhao/datasets/bignn/rand/rand1m3d_gold_rnn10.txt";
    const char *knnsrc3  = "/home/wlzhao/datasets/bignn/rand/rand1m3d_k=24v9.txt";
    const char *knngrd4  = "/home/wlzhao/datasets/bignn/rand/rand1m32d_gold_knn30.txt";
    const char *knnsrc4  = "/home/wlzhao/datasets/bignn/rand/rand1m32d_k=24v10.txt";
    const char *knngrd5  = "/home/wlzhao/datasets/bignn/rand/rand1m100d_gold_knn40.txt";
    const char *knnsrc5  = "/home/wlzhao/datasets/bignn/rand/rand1m100d_k=40v14.txt";
    const char *knngrd6  = "/home/wlzhao/datasets/bignn/rand/rand1m8d_gold_knn20x.txt";
    const char *knnsrc6  = "/home/wlzhao/datasets/bignn/rand/rand1m8d_k=8v10.txt";
    const char *knngrd7  = "/home/wlzhao/datasets/bignn/glove/glove1m_gold_knn40.txt";
    const char *knnsrc7  = "/home/wlzhao/datasets/bignn/glove/glove1m_k=40v13.txt";
    const char *knngrd8  = "/home/wlzhao/datasets/bignn/nuswide/nusw_gold_kiknn40.txt";
    const char *knnsrc8  = "/home/wlzhao/datasets/bignn/nuswide/nusw_kik=40v13.txt";
    const char *knngrd9  = "/home/wlzhao/datasets/bignn/rand/rand100k2d_gold_knn30.txt";
    const char *knnsrc9  = "/home/wlzhao/datasets/bignn/rand/rand100k2d_k=6v9.txt";
    const char *knngrd10  = "/home/wlzhao/datasets/bignn/rand/rand100k100d_gold_knn30.txt";
    const char *knnsrc10  = "/home/wlzhao/datasets/bignn/rand/rand100k100d_k=42v13.txt";
    const char *knngrd11  = "/home/wlzhao/datasets/bignn/rand/rand10m3d_gold_knn30.txt";
    const char *knnsrc11  = "/home/wlzhao/datasets/bignn/rand/rand10m3d_k=30v5.txt";
    const char *knngrd12  = "/home/wlzhao/datasets/bignn/rand/rand100k2d_gold_knn30.txt";
    const char *knnsrc12  = "/home/wlzhao/datasets/bignn/rand/rand100k2d_k=42v13.txt";
    const char *knngrd13  = "/home/wlzhao/datasets/bignn/yfcc1m/yfcc1m_gold_knn30.txt";
    const char *knnsrc13  = "/home/wlzhao/datasets/bignn/yfcc1m/yfcc1m_k=40v9.txt";
    const char *knngrd14  = "/home/wlzhao/datasets/bignn/sift10m/sift10m_gold_knn40.txt";
    const char *knnsrc14  = "/home/wlzhao/datasets/bignn/sift10m/sift10m_k=40v13.txt";
    const char *knngrd15  = "/home/wlzhao/datasets/odd/s1feats1m_ref_gold_knn40.txt";
    const char *knnsrc15  = "/home/wlzhao/datasets/odd/s1feats1m_ref_k=40v10.txt";
    const char *knngrd16  = "/home/wlzhao/datasets/bignn/rand/grand100k50d_gold_knn40.txt";
    const char *knnsrc16  = "/home/wlzhao/datasets/bignn/rand/grand100k50d_kgraph_k=10.txt";

      const char *knnsrc17 = "/home/pclin/osptag/graph_out.txt";
      ///const char *knnsrc18 = "/home/pclin/nns_benchmark-master/distributed/rand/test/sift1m_40.txt";

    const char *knngrd17 = "/home/pclin/osptag/rand10m8d_K20_gold.txt";

    ///const char *knngrd18 = "/home/pclin/nns_benchmark-master/data/sift/sift1m_gold_knn40.txt";

    const char *result1  = "/home/wlzhao/datasets/bignn/rslt/sift1m_rq8v10.txt";
    const char *result2  = "/home/wlzhao/datasets/bignn/rslt/sift1m_knn_rq8.txt";
    const char *result3  = "/home/wlzhao/datasets/bignn/rslt/efanna_rslt.txt";
    const char *result4  = "/home/wlzhao/datasets/bignn/rslt/efanna_rslt.txt";
    const char *gistrslt = "/home/wlzhao/datasets/bignn/rslt/qnn_gist_rslt.txt";
    const char *siftrslt = "/home/wlzhao/datasets/bignn/rslt/qnn_sift_rslt.txt";
    const char *kgraphrslt = "/home/wlzhao/datasets/bignn/rslt/sift1m_kgraph_rslt.txt";
    const char *nuswrslt = "/home/wlzhao/datasets/bignn/rslt/nusw_pq50.txt";
    const char *randrslt = "/home/wlzhao/datasets/bignn/rslt/flann_result_top1000.txt";
    const char *gloverslt = "/home/wlzhao/datasets/bignn/rslt/qnn_glove_rslt.txt";
    const char *sift100rslt = "/home/wlzhao/datasets/bignn/rslt/qnn_sift100k_rslt.txt";
    const char *yfccrslt = "/home/wlzhao/datasets/bignn/rslt/qnn_yfcc1m_rslt.txt";
    const unsigned int topk1 = 8192;
    const unsigned int topk2 = 10;
    const unsigned int tops[12] = {1, 2, 4, 5, 8, 10, 16, 32, 50, 64, 100, 128};

    ///Evaluator::getRGraphRecall(knnsrc3, knngrd3);
    ///Evaluator::getKGraphRecall(knnsrc1, knngrd1, 1, topk2);
    ///Evaluator::getKGraphRecall(knnsrc2, knngrd2, 1, topk2);
    ///Evaluator::getKGraphRecall(knnsrc3, knngrd3, 1, topk2);
    ///Evaluator::getKGraphRecall(knnsrc4, knngrd4, 1, topk2);
    ///Evaluator::getKGraphRecall(knnsrc5, knngrd5, 1, topk2);
    ///Evaluator::getKGraphRecall(knnsrc6, knngrd6, 1, topk2);
    ///Evaluator::getKGraphRecall(knnsrc7, knngrd7, 1, topk2);
    ///Evaluator::getKGraphRecall(knnsrc8, knngrd8, 1, topk2);
    ///Evaluator::getKGraphRecall(knnsrc9, knngrd9, 1, topk2);
    ///Evaluator::getKGraphRecall(knnsrc18, knngrd18, 10, 1);
    ///Evaluator::getKGraphRecall(knnsrc18, knngrd18, 10, 10);
    Evaluator::getKGraphRecall(knnsrc17, knngrd17, 1);
    Evaluator::getKGraphRecall(knnsrc17, knngrd17, 10);
    ///Evaluator::getKGraphRecall(knnsrc11, knngrd11, 1, topk2);
    ///Evaluator::getKGraphRecall(knnsrc12, knngrd12, 1, topk2);
    ///Evaluator::getKGraphRecall(knnsrc13, knngrd13, 1, topk2);
    ///Evaluator::getKGraphRecall(knnsrc14, knngrd14, 1, topk2);
    ///Evaluator::getKGraphRecall(knnsrc15, knngrd15, 1, topk2);
    ///Evaluator::getKGraphRecall(knnsrc16, knngrd16, 1, topk2);

    ///Evaluator::NNRecall(gistrlst, gistgrd, cout, 1);
    ///Evaluator::NNRecall(result1, siftgrd, cout, 1);
    //for(int i = 4; i < 5; i++)
    {
        ///Evaluator::NNRecall(gistrslt, gistgrd, cout, tops[5]);
        ///Evaluator::NNRecall(kgraphrslt, siftgrd, cout, 10);
        ///Evaluator::NNRecall(gloverslt, glovegrd, cout, 10);
        ///Evaluator::NNRecall(nuswrslt, nuswgrd, cout, 16);
        ///Evaluator::NNRecall(yfccrslt, yfccgrd, cout, 10);
        ///Evaluator::NNRecall(randrslt, randgrd, cout, 10);
        ///Evaluator::NNRecall(sift100rslt, sift100kgrd, cout, 10);

        //Evaluator::NNRecall(result4, gistgrd, cout, tops[i]);
    }

    return ;
}

