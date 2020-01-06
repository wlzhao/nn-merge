#include <iostream>
#include <cstring>
#include <cassert>
#include <iomanip>
#include <cstdio>
#include <cmath>

#include "iodelegator.h"
#include "vstring.h"
#include "vmath.h"
#include "distopt.h"

using namespace std;

long IODelegator::getFileSize(const char *srcfn)
{
    ifstream *inStrm = new ifstream(srcfn, ios::in|ios::binary);
    inStrm->seekg(0, ios::end);
    long sz = (long)inStrm->tellg();
    inStrm->close();
    delete inStrm;

    return sz;
}

float *IODelegator::loadTruncMat(const char *srcFn, unsigned int &row, unsigned int &dim, vector<unsigned int> &ids)
{
    ifstream *inStrm = new ifstream(srcFn, ios::in);
    unsigned int i = 0, j = 0, id;
    (*inStrm)>>row;
    (*inStrm)>>dim;
    float *mat = new float[dim*row];
    for(i = 0; i < row; i++)
    {
         (*inStrm)>>id;
         ids.push_back(id);
         for(j = 0; j < dim; j++)
         {
              (*inStrm)>>mat[i*dim+j];
         }
    }
    inStrm->close();
    return mat;
}

void IODelegator::saveTruncMat(const char *dstFn, const unsigned int dim, float *mat, set<unsigned int> &ids)
{
    ofstream *outStrm = new ofstream(dstFn, ios::out);
    set<unsigned int>::iterator sit;
    unsigned int row = ids.size(), i = 0, j = 0;
    (*outStrm)<<row<<" "<<dim<<endl;
    for(sit = ids.begin(), i = 0; sit != ids.end(); sit++,i++)
    {
         (*outStrm)<<*sit;
         for(j = 0; j < dim; j++)
         {
              (*outStrm)<<" "<<mat[i*dim+j];
         }
         (*outStrm)<<endl;
    }
    outStrm->close();
}


unsigned int *IODelegator::loadKNNGraph(const char *srcFn, unsigned int &gRow, unsigned int &gDim, const unsigned int topk0)
{
    unsigned int *knnGraph = NULL, r = 0, i = 0, idx = 0, tmp = 0, k;
    ifstream *inStrm = new ifstream(srcFn, ios::in);
    if(!inStrm->is_open())
    {
        cout<<"File '"<<srcFn<<"' cannot open for read!\n";
        exit(0);
    }
    gDim = 0;
    (*inStrm)>>gRow;
    (*inStrm)>>gDim;
    assert(gDim >= topk0);
    knnGraph = (unsigned int*)calloc(gRow, topk0*sizeof(unsigned int));
    while(!inStrm->eof() && r < gRow)
    {
        (*inStrm)>>tmp;
        (*inStrm)>>k;
        for(i = 0; i < k; i++)
        {
            if(i < topk0)
            {
                (*inStrm)>>knnGraph[idx];
                idx++;
            }
            else
            {
                (*inStrm)>>tmp;
            }
        }
        r++;
        ///cout<<"\r\r\r\r\t"<<r;
    }
    gDim = topk0;

    inStrm->close();
    return knnGraph;
}

unsigned int **IODelegator::loadANNGraph(const char *srcFn, unsigned int &gRow, unsigned int &gDim)
{
    unsigned int **annGraph = NULL;
    unsigned int n = 0, i = 0, id = 0, rid = 0;
    ifstream *inStrm = new ifstream(srcFn, ios::in);

    if(!inStrm->is_open())
    {
        cout<<"Open file '"<<srcFn<<"' failed!\n";
        return annGraph;
    }

    (*inStrm)>>gRow;
    (*inStrm)>>gDim;
    assert(gRow > 0);
    annGraph = new unsigned int*[gRow];
    //cout<<gRow<<"\t"<<dim<<endl;
    //memset(annGraph, 0, sizeof(unsigned int*)*gRow);
    while(!inStrm->eof())
    {
       (*inStrm)>>id;
       (*inStrm)>>n; n += 1;
       annGraph[id] = new unsigned int[n];
       annGraph[id][0] = n-1;
       for(i = 1; i < n; i++)
       {
           (*inStrm)>>rid;
           annGraph[id][i] = rid;
       }
    }
    inStrm->close();

    return annGraph;
}

void   IODelegator::saveTruncKnnGraph(const char *knnFn, const char *dstFn, unsigned int *idx, unsigned int idxNum)
{
    ofstream *outStrm = new ofstream(dstFn, ios::out);
    unsigned int i, j;

    unsigned int kRow, kDim, *knnGraph, crntIdx;
    knnGraph = IODelegator::loadKNNGraph(knnFn, kRow, kDim, 25);
    for(i = 0; i < idxNum; i++)
    {
        crntIdx = i*1000000;
        *outStrm << crntIdx << " ";
        for(j = 0; j < kDim; j++)
        {
            *outStrm << knnGraph[crntIdx*kDim + j] << "  ";
        }
        *outStrm << endl;
    }
    outStrm->close();
}

float *IODelegator::load_refSet(const char *srcFn, unsigned int &d, unsigned int &r)
{
    float *refMat = NULL;

    vector<unsigned int> ids;

    if(VString::endWith(srcFn, ".fvecs"))
    {
        refMat = IODelegator::load_fvecs(srcFn, d, r);
    }
    else if(VString::endWith(srcFn, ".bvecs"))
    {
        refMat = IODelegator::loadBvec10m(srcFn, d, r);
    }else if(VString::endWith(srcFn, ".tvecs"))
    {
        refMat = IODelegator::loadTruncMat(srcFn, r, d, ids);
    }
    else if(VString::endWith(srcFn, ".txt"))
    {
        refMat = IODelegator::loadMatrix(srcFn, r, d);
    }
    ids.clear();
    return refMat;
}



float  *IODelegator::loadMatrix(const char *fn, unsigned int &row, unsigned int &col)
{
    unsigned int irow = 0, idim = 0, loc = 0;
    float vals[2] = {0};
    assert(fn);

    ifstream *inStrm = new ifstream();
    inStrm->open(fn, ios::in);
    if(inStrm->fail())
    {
        cout<<"Fail to read "<<fn<<endl;
        delete inStrm;
        exit(0);
    }

    (*inStrm)>>vals[0];
    (*inStrm)>>vals[1];
    ///cout<<vals[0]<<"\t"<<vals[1]<<endl;
    row = (int)round(vals[0]);
    col = (int)round(vals[1]);
    float *mat = new float[row*col];


    for(irow = 0; irow < row; irow++)
    {
        loc = irow*col;
        for(idim = 0; idim < col; idim++)
        {
            (*inStrm) >>mat[loc + idim];
        }
    }
    inStrm->close();
    delete inStrm;
    return mat;
}

float *IODelegator::load_fvecs(const char *fvecFn, unsigned int &d, unsigned int &r)
{
    float *mat = NULL, *vect = NULL, *ppmat = NULL;
    unsigned bfsz = 0, bg = 0;
    unsigned int di = 0, line = 0;

    ifstream *inStrm = new ifstream(fvecFn, ios::in|ios::binary);
    if(!inStrm->is_open())
    {
        cout<<"File '"<<fvecFn<<"' cannot open for read!\n";
        return NULL;
    }
    r = d = 0;
    bg = inStrm->tellg();
    inStrm->read((char*)&di,  sizeof(unsigned int));
    d = di;
    bfsz = d*sizeof(float);
    inStrm->seekg(0, ios::end);
    r = ((long)inStrm->tellg() - bg)/(bfsz + sizeof(unsigned int));
    inStrm->close();

    mat   = new float[r*d];
    vect  = new float[d];
    ppmat = mat;

    inStrm = new ifstream(fvecFn, ios::in|ios::binary);
    while(!inStrm->eof() && line < r)
    {
        inStrm->read((char*)&di, sizeof(int));
        assert(di == d);

        bfsz = d*sizeof(float);
        inStrm->read((char*)vect, bfsz);

        line++;
        memcpy(ppmat, vect, bfsz);
        ppmat = ppmat + d;
    }

    delete [] vect;
    inStrm->close();

    return mat;
}

float *IODelegator::load_bvecs(const char *bvecFn, unsigned int &d, unsigned int &r)
{
    float *mat = NULL, *fvect = NULL, *ppmat = NULL;
    unsigned int di = 0, i = 0;
    unsigned long bfsz = 0, bg = 0;
    unsigned char *vect = NULL;

    ifstream *inStrm = new ifstream(bvecFn, ios::in|ios::binary);
    if(!inStrm->is_open())
    {
        cout<<"File '"<<bvecFn<<"' cannot open for read!\n";
        return NULL;
    }
    r = d = 0;
    bg = inStrm->tellg();
    inStrm->read((char*)&di,  sizeof(unsigned int));
    d = di;
    bfsz = d*sizeof(unsigned char);
    inStrm->seekg(0, ios::end);
    r = ((long)inStrm->tellg() - bg)/(bfsz + sizeof(unsigned int));
    inStrm->close();

    mat   = new float[r*d];
    vect  = new unsigned char[d];
    fvect = new float[d];
    ppmat = mat;

    inStrm = new ifstream(bvecFn, ios::in|ios::binary);
    while(!inStrm->eof())
    {
        inStrm->read((char*)&di, sizeof(int));
        assert(di == d);

        bfsz = d*sizeof(unsigned char);
        inStrm->read((char*)vect, bfsz);
        for(i = 0; i < d; i++)
        {
            fvect[i] = (unsigned char)vect[i];
        }

        memcpy(ppmat, fvect, d*sizeof(float));
        ppmat = ppmat + d;
    }

    delete [] vect;
    delete [] fvect;
    fvect = NULL;
    vect  = NULL;
    inStrm->close();
    return mat;
}

float *IODelegator::loadBvec10m(const char *bvecfn, unsigned int &d, unsigned int &r)
{
    unsigned char *vect = NULL;
    float *mat = NULL, *fvect = NULL, *ppmat = NULL;
    unsigned int di = 0, ir = 0;
    d = 0;
    r = 0;
    unsigned long bg = 0, fsize = 0, bfsz = 0;

    ifstream *inStrm = new ifstream(bvecfn, ios::in|ios::binary);
    if(!inStrm->is_open())
    {
        cout<<"File '"<<bvecfn<<"' cannot open for read!\n";
        return NULL;
    }

    //inStrm = new ifstream(fvecfn, ios::in|ios::binary);
    bg = (long)inStrm->tellg();
    inStrm->read((char*)&di,  sizeof(unsigned int));
    d = di;
    bfsz = d*sizeof(unsigned char);
    inStrm->seekg(0, ios::end);
    fsize = (unsigned long)inStrm->tellg() - bg;
    inStrm->close();
    r = fsize/(sizeof(unsigned int) + bfsz);
    if(r == 0)
    {
        cout<<"No data has been loaded!\n";
        r = d = 0;
        return NULL;
    }
    r = r < 10000000? r:10000000;
    vect   = new unsigned char[d];
    mat    = new float[r*d];
    fvect  = new float[d];
    memset(mat, 0, sizeof(float)*r*d);
    ppmat  = mat;
    inStrm = new ifstream(bvecfn, ios::in|ios::binary);
    int i;
    while(!inStrm->eof() && ir < r)
    {
        inStrm->read((char*)&di,  sizeof(unsigned int));

        if(di == 0)
            continue;

        bfsz = di*sizeof(unsigned char);
        inStrm->read((char*)vect, bfsz);
        for(i = 0; i < di; i++)
        {
            fvect[i] = (unsigned char)vect[i];
            //cout << (int)vect[i] << " " << fvect[i] << " " << endl;
        }
        memcpy(ppmat, fvect, di*sizeof(float));
        memset(vect, 0, bfsz);
        ppmat = ppmat + d;
        di = 0;
        ir++;
    }

    delete [] vect;
    delete [] fvect;
    vect = NULL;
    fvect = NULL;
    inStrm->close();
    return mat;
}

float *IODelegator::loadBvec100m(const char *bvecfn, unsigned int &d, unsigned int &r)
{
    unsigned char *vect = NULL;
    float *mat = NULL, *fvect = NULL, *ppmat = NULL;
    unsigned int di = 0, ir = 0;
    d = 0;
    r = 0;
    unsigned long bg = 0, fsize = 0, bfsz = 0;

    ifstream *inStrm = new ifstream(bvecfn, ios::in|ios::binary);
    if(!inStrm->is_open())
    {
        cout<<"File '"<<bvecfn<<"' cannot open for read!\n";
        return NULL;
    }

    //inStrm = new ifstream(fvecfn, ios::in|ios::binary);
    bg = (long)inStrm->tellg();
    inStrm->read((char*)&di,  sizeof(unsigned int));
    d = di;
    bfsz = d*sizeof(unsigned char);
    inStrm->seekg(0, ios::end);
    fsize = (unsigned long)inStrm->tellg() - bg;
    inStrm->close();
    r = fsize/(sizeof(unsigned int) + bfsz);
    r = r > 100000000? 100000000:r;
    if(r == 0)
    {
        cout<<"No data has been loaded!\n";
        r = d = 0;
        return NULL;
    }
    vect   = new unsigned char[d];
    mat    = (float*)calloc(r, sizeof(float)*d);
    fvect  = new float[d];
    ppmat  = mat;
    inStrm = new ifstream(bvecfn, ios::in|ios::binary);
    int i;
    while(!inStrm->eof() && ir < r)
    {
        inStrm->read((char*)&di,  sizeof(unsigned int));

        if(di == 0)
            continue;

        bfsz = di*sizeof(unsigned char);
        inStrm->read((char*)vect, bfsz);
        for(i = 0; i < di; i++)
        {
            fvect[i] = (unsigned char)vect[i];
            //cout << (int)vect[i] << " " << fvect[i] << " " << endl;
        }
        memcpy(ppmat, fvect, di*sizeof(float));
        memset(vect, 0, bfsz);
        ppmat = ppmat + d;
        di = 0;
        ir++;
    }

    delete [] vect;
    delete [] fvect;
    vect = NULL;
    fvect = NULL;
    inStrm->close();
    return mat;
}


void IODelegator::loadSMat(const char *srcFn, unsigned int &r, vector<unordered_set<unsigned> > &smat)
{
    ifstream inStrm;
    unsigned int nRow = 0, num = 0, i = 0, ir = 0;
    unsigned int val = 0;
    inStrm.open(srcFn, ios::in);

    if(!inStrm.is_open())
    {
        cout<<"File '"<<srcFn<<"' cannot open for read!\n";
        r = 0;
        return ;
    }
    inStrm>>nRow;
    smat.resize(nRow);
    r = nRow;
    while(!inStrm.eof() && ir < nRow)
    {
        inStrm>>num;
        unordered_set<unsigned> &onerow = smat[ir];
        for(i = 0; i < num; i++)
        {
            inStrm>>val;
            onerow.insert(val);
        }
        ir++;
    }

    inStrm.close();
}

float *IODelegator::bvec2txt(const char *bvecfn, const char *dstfn, unsigned int &d, unsigned int &r)
{
    unsigned char *vect = NULL;
    float *mat = NULL, *fvect = NULL, *ppmat = NULL;
    unsigned int di = 0, ir = 0;
    d = 0;
    r = 0;
    unsigned long bg = 0, fsize = 0, bfsz = 0;
    ofstream dstOstrm(dstfn, ios::out);
    {
        if(!dstOstrm.is_open())
        {
            cout << dstfn << " can't open" << endl;
            exit(0);
        }
    }

    ifstream *inStrm = new ifstream(bvecfn, ios::in|ios::binary);
    if(!inStrm->is_open())
    {
        cout<<"File '"<<bvecfn<<"' cannot open for read!\n";
        return NULL;
    }

    //inStrm = new ifstream(fvecfn, ios::in|ios::binary);
    bg = (long)inStrm->tellg();
    inStrm->read((char*)&di,  sizeof(unsigned int));
    d = di;
    bfsz = d*sizeof(unsigned char);
    inStrm->seekg(0, ios::end);
    fsize = (unsigned long)inStrm->tellg() - bg;
    inStrm->close();
    r = fsize/(sizeof(unsigned int) + bfsz);
    if(r == 0)
    {
        cout<<"No data has been loaded!\n";
        r = d = 0;
        return NULL;
    }
    vect   = new unsigned char[d];
    mat    = new float[r*d];
    fvect  = new float[d];
    memset(mat, 0, sizeof(float)*r*d);
    ppmat  = mat;
    inStrm = new ifstream(bvecfn, ios::in|ios::binary);
    int i;
    dstOstrm << r << "  " << d << endl;
    //cout << "sizeof unsign char " << sizeof(int) << endl;
    while(!inStrm->eof() && ir < r)
    {
        inStrm->read((char*)&di,  sizeof(unsigned int));

        if(di == 0)
            continue;

        bfsz = di*sizeof(unsigned char);
        inStrm->read((char*)vect, bfsz);
        for(i = 0; i < di; i++)
        {
            fvect[i] = (int)vect[i];
            dstOstrm << fvect[i] << " ";
            //cout << (int)vect[i] << " " << fvect[i] << " " << endl;
        }
        dstOstrm << endl;
        memcpy(ppmat, fvect, di*sizeof(float));
        memset(vect, 0, bfsz);
        ppmat = ppmat + d;
        di = 0;
        ir++;
    }

    dstOstrm.close();
    delete [] vect;
    delete [] fvect;
    vect = NULL;
    fvect = NULL;
    inStrm->close();
    return mat;
}

float *IODelegator::loadXQVocab(const char *vocabFn,  unsigned int &nqDim, unsigned int numb[16], unsigned int &dv0)
{
    float *mat = NULL;
    unsigned int i = 0, k = 0, idx = 0, sum = 0;
    float tmp = 0;

    ifstream *inStrm = new ifstream(vocabFn, ios::in);
    if(!inStrm->is_open())
    {
        cout<<"File "<<vocabFn<<" cannot open for read!"<<endl;
        exit(0);
    }

    (*inStrm)>>dv0;
    for(i = 0; i < dv0; i++)
    {
        (*inStrm)>>numb[i];
        sum += numb[i];
        ///cout<<numb[i]<<" ";
    }
    ///cout<<endl;

    (*inStrm)>>nqDim;

    mat = new float[sum*nqDim];

    for(i = 0; i < sum; i++)
    {
        for(k = 0; k < nqDim; k++, idx++)
        {
            (*inStrm)>>tmp;
            mat[idx] = tmp;
        }
    }

    inStrm->close();

    return mat;
}

float *IODelegator::loadVocab(const char *fn, unsigned int &row, unsigned int &col)
{
    assert(fn);
    float vals[3] = {0};

    ifstream *inStrm = new ifstream();
    inStrm->open(fn, ios::in);
    if(inStrm->fail())
    {
        cout<<"Fail to read "<<fn<<endl;
        delete inStrm;
        exit(0);
    }

    (*inStrm)>>vals[0];
    (*inStrm)>>vals[1];
    (*inStrm)>>vals[2];
    row = (int)round(vals[1]);
    col = (int)round(vals[2]);
    float *mat = new float[row*col];

    unsigned int irow = 0, idim = 0, loc = 0;

    for(irow = 0; irow < row; irow++)
    {
        loc = irow*col;
        for(idim = 0; idim < col; idim++)
        {
            (*inStrm) >>mat[loc + idim];
        }
    }
    inStrm->close();
    delete inStrm;
    return mat;
}

void IODelegator::loadMultiMap(const char *srcfn, vector<vector<unsigned int> > &sup_map)
{
    unsigned int val, size, numb, j = 0;
    ifstream *inStrm = new ifstream(srcfn);

    (*inStrm)>>numb;
    for(j = 0; j < numb; j++)
    {
        vector<unsigned int> row;
        sup_map.push_back(row);
    }

    while(!inStrm->eof())
    {
        (*inStrm)>>size;
        (*inStrm)>>j;
        vector<unsigned> &crntVect = sup_map[j];

        for(unsigned int i = 0; i < size; i++)
        {
            (*inStrm)>>val;
            crntVect.push_back(val);
        }
    }
    inStrm->close();

    return ;
}

bool IODelegator::loadWghs(const char *weightfn, map<unsigned int,float> &i2wMap,
                           const unsigned int offset)
{
    assert(weightfn);
    FILE *fp = fopen(weightfn,"r");
    if(fp== NULL)
    {
        cout<<"Open file "<<weightfn<<" error!\n";
        exit(1);
    }

    char oneline[LONGCHAR];
    unsigned int IMG_ID;
    int num = 0, nargs = 0;
    float weight = 0;
    ///cout<<weightfn<<endl;
    char *ptStr = NULL;

    while(!feof(fp))
    {
        ptStr = fgets(oneline, LONGCHAR, fp);

        if(ptStr == NULL)
        continue;

        if(!strcmp(ptStr, ""))
		continue;

        if(strcmp(oneline, ""))
        {
            nargs = sscanf(oneline,"%d %f %d",&IMG_ID, &weight, &num);

            if(nargs != 3)
            continue;

            IMG_ID = IMG_ID + offset;
            i2wMap.insert(pair<unsigned int, float>(IMG_ID, weight));
            strcpy(oneline,"");
        }
    }

    fclose(fp);

    return true;
}

bool IODelegator::loadWghs(const char *weightfn, map<unsigned int,float> &i2wMap, map<unsigned int, unsigned int> &pmap,
                           map<unsigned int, unsigned int>  &_pmap, map<unsigned int, int>  &tmap, const unsigned int kf_offset)
{
    assert(weightfn);
    FILE *fp = fopen(weightfn,"r");
    if(fp== NULL)
    {
        cout<<"Open file "<<weightfn<<" error!\n";
        exit(1);
    }

    ///cout<<"i am here\n";

    char oneline[LONGCHAR];
    unsigned int IMG_ID, crnt_vid, pre_vid = 0, tm_code;
    int nargs = 0;
    float weight = 0;
    unsigned vid_offset =  pmap.size();
    char *ptStr = NULL;

    while(!feof(fp))
    {
        ptStr = fgets(oneline, LONGCHAR, fp);

        if(ptStr == NULL)
        continue;

        if(!strcmp(ptStr, ""))
        continue;

        if(strcmp(oneline, ""))
        {
            nargs = sscanf(oneline,"%d %d %d %f",&IMG_ID, &crnt_vid, &tm_code, &weight);

            if(nargs != 4)
            continue;

            IMG_ID = IMG_ID + kf_offset;
            _pmap.insert(pair<unsigned int,unsigned int>(IMG_ID, crnt_vid+vid_offset));
            if(crnt_vid != pre_vid)
            {
                pmap.insert(pair<unsigned int,unsigned int>(crnt_vid+vid_offset, IMG_ID));
            }
            i2wMap.insert(pair<unsigned int,float>(IMG_ID,weight));
            tmap.insert(pair<unsigned int,int>(IMG_ID,tm_code));
            strcpy(oneline,"");
            pre_vid = crnt_vid;
        }
    }

    fclose(fp);

    return true;
}

void IODelegator::extrTrucMat(const unsigned int numb, const char *srcFn, const char *dstFn)
{
   map<unsigned int, unsigned int> idx2key;
   vector<unsigned int> ids;
   set<unsigned int> seeds;

   unsigned int i = 0, dim = 0, nRow = 0, loc = 0;
   unsigned int counts = 0;
   float *mat = NULL;

   cout<<"Load matrix ........................ ";
   if(VString::endWith(srcFn, ".txt"))
   {
       mat = IODelegator::loadMatrix(srcFn, nRow, dim);
       for(i = 0; i < nRow; i++)
       {
          ids.push_back(i);
       }
   }else if(VString::endWith(srcFn, ".tvecs"))
   {
       mat = IODelegator::loadTruncMat(srcFn, nRow, dim, ids);
   }else if(VString::endWith(srcFn, ".fvecs"))
   {
       mat = IODelegator::load_fvecs(srcFn, dim, nRow);
       for(i = 0; i < nRow; i++)
       {
          ids.push_back(i);
       }
   }else if(VString::endWith(srcFn, ".bvecs"))
   {
       mat = IODelegator::loadBvec10m(srcFn, dim, nRow);
       for(i = 0; i < nRow; i++)
       {
          ids.push_back(i);
       }
   }

   cout<<nRow<<"\t"<<dim<<endl;
   VMath::randSeed(seeds, nRow-1, numb);

   cout<<"Extracting seed vectors ........... ";
   ofstream *outStrm = new ofstream(dstFn, ios::out);
   (*outStrm)<<seeds.size()<<"\t"<<dim<<endl;
   for(auto sit = seeds.begin(); sit != seeds.end(); sit++)
   {
        loc = *sit*dim;
        (*outStrm)<<ids[*sit];
        ///cout<<*sit<<"\t"<<ids[*sit]<<endl;
        assert(*sit < nRow);
        for(i = 0; i < dim; i++)
        {
            (*outStrm)<<" "<<mat[loc+i];
        }
        (*outStrm)<<endl;
        counts++;
   }
   outStrm->close();
   cout<<endl;
   delete [] mat;
   mat = NULL;
   cout<<"Done: "<<counts<<endl;
   seeds.clear();
   ids.clear();
}

void IODelegator::extrTrucMat(const char *srcTrc, const char *srcFn, const char *dstFn)
{
   set<unsigned int> seeds;
   unsigned int i = 0, j = 0, dim, nRow = 0, loc = 0;

   ifstream *inStrm = new ifstream(srcTrc, ios::in);
   if(!inStrm->is_open())
   {
      cout<<"File "<<srcTrc<<" cannot open for writting!\n";
      exit(0);
   }
   while(!inStrm->eof() && nRow < 10000)
   {
       (*inStrm)>>i;
       (*inStrm)>>j;
       seeds.insert(j);
       nRow++;
   }
   //cout<<i<<"\t"<<j<<endl;
   inStrm->close();
   float *mat = IODelegator::loadMatrix(srcFn, nRow, dim);
   int counts = 0;
   cout<<nRow<<"\t"<<dim<<endl;
   ofstream *outStrm = new ofstream(dstFn, ios::out);
   set<unsigned int>::iterator sit;
   for(sit = seeds.begin(); sit != seeds.end(); sit++)
   {
        loc = *sit*dim;
        (*outStrm)<<*sit;
        for(i = 0; i < dim; i++)
        {
            (*outStrm)<<" "<<mat[loc+i];
        }
        (*outStrm)<<endl;
        counts++;
   }
   outStrm->close();
   cout<<endl;
   delete [] mat;
   mat = NULL;
   cout<<"OK: "<<counts<<endl;
   seeds.clear();

}

void IODelegator::load_k2imap(const char *srcfn, map<string, unsigned int> &referTab)
{
    FILE *fp = fopen(srcfn, "r");

    if(fp == NULL)
        return ;

    char key[256];
    char *crnt_key = NULL;
    int klen = 0, label = 0, nargs = 0;

    while(!feof(fp))
    {
        nargs = fscanf(fp, "%s %d", key, &label);

        if(nargs != 2)
        continue;

        klen = strlen(key)+1;
        crnt_key = new char[klen];

        strcpy(crnt_key, key);

        referTab.insert(pair<const char*, unsigned int>(crnt_key, label));
    }
    fclose(fp);
    return ;
}

void IODelegator::load_k2imap(const char *srcfn, map<string, int> &referTab)
{
    FILE *fp = fopen(srcfn, "r");

    if(fp == NULL)
        return ;

    char key[256];
    char *crnt_key;
    int klen = 0, label = 0, nargs = 0;

    while(!feof(fp))
    {
        nargs = fscanf(fp, "%s %d", key, &label);

        if(nargs != 2)
        continue;

            klen = strlen(key) + 1;
        crnt_key = new char[klen];
        strcpy(crnt_key, key);

        referTab.insert(pair<const char*, int>(crnt_key, label));
    }
    fclose(fp);

    return ;
}

void IODelegator::load_i2kmap(const char *nmapfn, map<unsigned int, const char*> &i2kMap,
                              const unsigned int offset)
{
     assert(nmapfn);
    FILE *fp = fopen(nmapfn, "r");

    if(fp == NULL)
    {
        cout<<"Open file '"<<nmapfn<<"' error!\n";
        exit(1);
    }

    char oneline[LONGCHAR];
    char fname[FNLEN];

    unsigned int KF_ID, ID = 0, localID = 0;
    int nargs = 0;
    char *name = NULL;

    while(!feof(fp))
    {
        char *strPt = fgets(oneline, LONGCHAR, fp);

        if(strPt == NULL)
            continue;

        if(!strcmp(oneline, "\n")||!strcmp(oneline, ""))
            continue;

        if(oneline[0] == '#')
            continue;

        nargs = sscanf(oneline,"%d %s\n", &KF_ID, fname);

        if(nargs != 2)
        continue;

        name = new char[strlen(fname)+1];
        strcpy(name, fname);
        name[strlen(fname)] = '\0';
        localID++;
        ID = localID +  offset;

        i2kMap.insert(pair<unsigned int, const char*>(ID, name));
        strcpy(oneline, "");

    }
    fclose(fp);
    return ;
}

vector<const char*> IODelegator::loadNameList(const char *kflstfn, const char *msg)
{
    assert(kflstfn);
    FILE *fp = fopen(kflstfn,"r");

    vector<const char*> namelst;

    if(fp == NULL)
    {
        cout<<"Open file "<<kflstfn<<" error!\n";
        exit(1);
    }

    cout<<msg;

    char oneline[LONGCHAR];
    char fname[FNLEN];

    unsigned int IMG_ID = 0, i = 0;
    char *name = NULL, *ptStr = NULL;
    int nargs = 0;

    while(!feof(fp))
    {
        ptStr = fgets(oneline, LONGCHAR, fp);

        if(ptStr == NULL)
        continue;

        if(!strcmp(ptStr, ""))
        continue;

        if(!strcmp(oneline, "\n")||!strcmp(oneline, ""))
        {
            continue;
        }

        if(oneline[0] != '#' || oneline[0] != '\n')
        {
            nargs = sscanf(oneline, "%d %s\n", &IMG_ID, fname);

            if(nargs != 2)
            continue;

             name = new char[strlen(fname)+1];
            strcpy(name,fname);

            name[strlen(fname)] = '\0';
            namelst.push_back(name);
            i++;
        }
        strcpy(oneline,"");
    }
    fclose(fp);
    cout<<i<<endl;

    return namelst;
}

vector<const char*> IODelegator::loadStrings(const char *kflstfn)
{
    assert(kflstfn);
    FILE *fp = fopen(kflstfn,"r");

    vector<const char*> namelst;

    if(fp == NULL)
    {
        cout<<"Open file "<<kflstfn<<" error!\n";
        exit(1);
    }

    cout<<"Loading Name List ......... ";

    char oneline[LONGCHAR];
    char fname[FNLEN];

    char *name = NULL, *ptStr = NULL;
    int i = 0;

    while(!feof(fp))
    {
        ptStr = fgets(oneline, LONGCHAR, fp);

        if(ptStr == NULL)
        continue;

        if(!strcmp(ptStr, ""))
        continue;

        if(!strcmp(oneline, "\n")||!strcmp(oneline, ""))
        {
            continue;
        }

        if(oneline[0] != '#' || oneline[0] != '\n')
        {
            sscanf(oneline, "%s\n", fname);
            name = new char[strlen(fname)+1];
            strcpy(name,fname);

            name[strlen(fname)] = '\0';
            namelst.push_back(name);
            i++;
        }
        strcpy(oneline, "");
    }
    fclose(fp);
    cout<<i<<endl;
    return namelst;
}

void IODelegator::load_map(const char *srcfn, map<string,unsigned char> &referTab)
{
    FILE *fp = fopen(srcfn, "r");
    if(fp == NULL)
        return ;

    char key[256];
    char *crnt_key = NULL;
    int  klen = 0, nargs = 0;

    while(!feof(fp))
    {
        nargs = fscanf(fp, "%s", key);

        if(nargs != 1)
        continue;

            klen = strlen(key)+1;
        crnt_key = new char[klen];

        strcpy(crnt_key,key);

        referTab.insert(make_pair(crnt_key,1));
    }
    fclose(fp);
}

unsigned int *IODelegator::loadTrucDPGraph(const char *dpgFn, unsigned int &row, unsigned int *&knnLocs,
                                           vector<unsigned int> &ids)
{
    ifstream is(dpgFn, ios::in);
    unsigned int i = 0, j = 0, idx = 0, numb, val = 0;
    unsigned long cnts = 0, k = 0;
    is>>row;
    knnLocs = new unsigned int[row];
    knnLocs[0] = 0;
    for (i = 0; i < row; ++i)
    {
        is>>idx;
        is>>numb;
        cnts = cnts + numb + 1;
        for (j = 0; j < numb; ++j)
        {
           is>>val;
        }
        if(i+1 < row)
        knnLocs[i+1] = knnLocs[i] + numb+1;
    }
    is.close();

    ifstream is1(dpgFn, ios::in);
    unsigned int *knnGraph = new unsigned int[cnts];
    is1>>row;

    for (i = k = 0; i < row; ++i)
    {
        is1>>idx;
        is1>>numb;
        knnGraph[k] = numb; k++;
        ids.push_back(idx);
        for(j = 0; j < numb; ++j, k++)
        {
           is1>>knnGraph[k];
        }
    }
    is1.close();

    return knnGraph;
}

unsigned int *IODelegator::loadDPGraph(const char *dpgFn, unsigned int &row, unsigned int *&knnLocs)
{
    ifstream is(dpgFn, ios::in);
    unsigned int i = 0, j = 0, idx = 0, numb, val = 0;
    unsigned long cnts = 0, k = 0;
    assert(is.is_open());
    is>>row;
    knnLocs = new unsigned int[row];
    ///cout<<row<<endl;
    knnLocs[0] = 0;
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
            knnLocs[i+1] = knnLocs[i] + numb+1;
        }
        idx = 0; numb = 0;
    }
    is.close();

    ifstream is1(dpgFn, ios::in);
    unsigned int *knnGraph = new unsigned int[cnts];
    is1>>row;

    for (i = k = 0; i < row&& !is.eof(); ++i)
    {
        is1>>idx;
        is1>>numb;
        knnGraph[k] = numb; k++;
        for(j = 0; j < numb; ++j, k++)
        {
           is1>>knnGraph[k];
        }
    }
    is1.close();

    return knnGraph;
}


void IODelegator::loadSeeds(const char *srcFn, vector<unsigned int> &mySeeds)
{
     mySeeds.clear();
     unsigned int id, val, cnt = 0;
     ifstream *inStrm = new ifstream(srcFn, ios::in);
     while(!inStrm->eof())
     {
         (*inStrm)>>id;
         (*inStrm)>>val;
         mySeeds.push_back(id);
         cnt++;
         if(cnt == 100000)
         break;

     }
     inStrm->close();
}

void IODelegator::load_word2ref(const char *srcFn, unsigned int w2fMap[128])
{
     ifstream *inStrm = new ifstream(srcFn, ios::in);
     unsigned int nRow, val, key;
     (*inStrm)>>nRow;
     for(unsigned int  i = 0; i < nRow; i++)
     {
         (*inStrm)>>key;
         (*inStrm)>>val;
         w2fMap[key] = val;
     }
     inStrm->close();
}

unsigned int *IODelegator::loadDPG(const char *dpgFn, unsigned int &row, unsigned int &dim)
{
    ifstream is(dpgFn, ios::binary);
    char magic[8];
    uint32_t major;
    uint32_t minor;
    uint32_t N;
    is.read(magic, sizeof(magic));
    is.read(reinterpret_cast<char *>(&major), sizeof(major));
    is.read(reinterpret_cast<char *>(&minor), sizeof(minor));
    is.read(reinterpret_cast<char *>(&N), sizeof(N));
    row = N;
    dim = 195;

    multimap<unsigned int, unsigned int> maxK;
    multimap<unsigned int, unsigned int>::iterator mit;
    unsigned int *knnGraph;
    //knnGraph = new unsigned int[N*136];
    knnGraph =  (unsigned int*)calloc(N, 195*sizeof(unsigned int));
    //out << N << endl;
    for (unsigned i = 0; i < N; ++i)
    {
        unsigned K;
        is.read(reinterpret_cast<char *>(&K), sizeof(K));
        //Kout << "K1 = " << K << " ";
        is.read(reinterpret_cast<char *>(&K), sizeof(K));
        //Kout << "K2 = " << K << endl;
        maxK.insert(pair<unsigned int, unsigned int>(K, i));
        knnGraph[i*195] = K;
        //out << knnGraph[i*195] << " ";
        for (unsigned j = 0; j < K; ++j)
        {
            is.read(reinterpret_cast<char *>(&knnGraph[i*195+j+1]), sizeof(unsigned int));
            //out << knnGraph[i*195+j+1] << " ";
        }
        //out << endl;
                //is.read(reinterpret_cast<char *>(&knn[0]), K * sizeof(knn[0]));
    }
    //out.close();

    /***
    for(mit = maxK.begin(); mit != maxK.end(); mit++)
    {
        Kout << mit->first << endl;
    }
    Kout.close();
    /***/
    return knnGraph;
}

void IODelegator::randSelect(vector<unsigned int> &pool, set<unsigned int> &sel, const unsigned int sz)
{
     vector<unsigned int>::iterator vit;
     vector<unsigned int> left;
     set<unsigned int>::iterator sit;

     unsigned int i = 0, s = 0;
     do{
         i = rand()%pool.size();
         s = pool[i];
         sel.insert(s);
     }while(sel.size()< sz);
     for(vit = pool.begin(); vit != pool.end(); vit++)
     {
         if(sel.find(*vit) == sel.end())
         {
            left.push_back(*vit);
         }
     }
     pool.clear();
     pool.insert(pool.begin(),left.begin(), left.end());
     left.clear();
}

void IODelegator::load_NNmap(const char* srcFn, map<unsigned int, unsigned int> &nnmap)
{
    ifstream inStrm;
    inStrm.open(srcFn, ios::in);
    unsigned int key = 0, val = 0;
    if(!inStrm.is_open())
    {
        cout<<"File '"<<srcFn<<"' cannot open for write!\n";
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

void IODelegator::procNUSWide(const char *srcFn, const char *qryFn, const char *refFn, const char *trnFn)
{
    unsigned int row = 0, col = 0, i = 0, j = 0, loc = 0;
    vector<unsigned int>::iterator vit;
    set<unsigned int>::iterator sit;
    vector<unsigned int> pool;
    set<unsigned int> qry;
    set<unsigned int> trn;

    float *mat = IODelegator::loadMatrix(srcFn, row, col);
    cout<<row<<"\t"<<col<<endl;
    for(i = 0; i < row; i++)
    {
        pool.push_back(i);
    }

    IODelegator::randSelect(pool, qry, 1000);
    IODelegator::randSelect(pool, trn, 2048);
    ofstream *outStrm = new ofstream(qryFn, ios::out);
    (*outStrm)<<qry.size()<<" "<<col<<endl;
    for(sit = qry.begin(); sit != qry.end(); sit++)
    {
        loc = *sit*col;
        for(j = 0; j < col; j++)
        {
            (*outStrm)<<mat[loc+j]<<" ";
        }
        (*outStrm)<<endl;
    }

    outStrm->close();

    outStrm = new ofstream(trnFn, ios::out);
    (*outStrm)<<trn.size()<<" "<<col<<endl;
    for(sit = trn.begin(); sit != trn.end(); sit++)
    {
        loc = *sit*col;
        for(j = 0; j < col; j++)
        {
            (*outStrm)<<mat[loc+j]<<" ";
        }
        (*outStrm)<<endl;
    }

    outStrm->close();

    outStrm = new ofstream(refFn, ios::out);
    (*outStrm)<<pool.size()<<" "<<col<<endl;
    for(vit = pool.begin(); vit != pool.end(); vit++)
    {
        loc = *vit*col;
        for(j = 0; j < col; j++)
        {
            (*outStrm)<<mat[loc+j]<<" ";
        }
        (*outStrm)<<endl;
    }

    outStrm->close();

    pool.clear();
    qry.clear();
    trn.clear();
}

void IODelegator::test()
{
    const char *srcFn  = "/home/wlzhao/datasets/bignn/nuswide/BoW_int.dat";
    const char *srcFn1 = "/home/wlzhao/datasets/bignn/glove/glove100k_rand.tvecs";
    const char *srcFn2 = "/home/wlzhao/datasets/bignn/gist/gist_base.fvecs";
    const char *srcFn3 = "/home/wlzhao/datasets/bignn/rand/rand1m100d.txt";
    const char *srcFn4 = "/home/wlzhao/datasets/bignn/annbench/kosarak27k_smat.txt";
    const char *qryFn  = "/home/wlzhao/datasets/bignn/nuswide/nusw_qry.txt";
    const char *refFn  = "/home/wlzhao/datasets/bignn/nuswide/nusw_26k.txt";
    const char *trnFn  = "/home/wlzhao/datasets/bignn/nuswide/nusw_trn.txt";
    const char *dstFn1 = "/home/wlzhao/datasets/bignn/glove/glove10k_rand.tvecs";
    const char *dstFn2 = "/home/wlzhao/datasets/bignn/gist/gist100k_rand.tvecs";
    const char *dstFn3 = "/home/wlzhao/datasets/bignn/rand/rand100k100d.trc";
    const char *dstFn4 = "/home/wlzhao/datasets/bignn/rand/rand1k100d.txt";
    const char *src = "/home/pclin/nns_benchmark-master/data/annbench/kosarak27k.smat";
    ///IODelegator::procNUSWide(srcFn, qryFn, refFn, trnFn);
    ///IODelegator::extrTrucMat(100000, srcFn3, dstFn3);
    vector<unordered_set<unsigned> > smat;
    unsigned int r = 0;
    IODelegator::loadSMat(src, r, smat);
    cout<<r<<"\t"<<smat.size()<<endl;
    cout<<"dist: "<<DistOpt::jaccard(smat[0], smat[3939])<<endl;
    cout<<"dist: "<<DistOpt::jaccard(smat[0], smat[53714])<<endl;
    cout<<"dist: "<<DistOpt::jaccard(smat[0], smat[74530])<<endl;
}


