#ifndef GVETEX_H
#define GVETEX_H
#include <vector>
struct GVetex {
public:
  float dst;        //4
  unsigned int id;  //4
  unsigned int rd;  //4
  unsigned int gd;
  unsigned char hits;  //1
  unsigned char xflag;
  GVetex *next; //8
  GVetex *prev; //8
  GVetex *rnxt; //8
  GVetex *rprv; //8

}__attribute__((__packed__)); ///for k-NN graph

struct HVetex {
public:
  unsigned int id;  //4
  unsigned int rd;  //4
  HVetex *next; //8
  HVetex *prev; //8
  HVetex *rnxt; //8
  HVetex *rprv; //8

}__attribute__((__packed__)); ///for sup k-NN graph

struct GHead {
public:
  GVetex *last;
  GVetex *head;
  unsigned int sz;
  unsigned int ht; //capacity
  float avgMMR;
  float mxDst;
}; ///for k-NN graph header

struct HHead {
public:
  HVetex *last;
  HVetex *head;
  unsigned int sz;
}; ///for sup k-NN graph header

struct linknode{
public:
    unsigned id;
    float dst;
    unsigned char flag;
    unsigned m;
    unsigned M;
    linknode *next;
    linknode *prev;
};

struct linklist{
public:
    linknode *head;
    linknode *last;
    unsigned sz;
    float max_dst;
};


struct RVetex{
public:
    float dst;         //4
    unsigned int id;   //4
    RVetex *next;      //8
}; ///for r-NN graph

struct RHead{
public:
  RVetex *head;
  unsigned int sz;
  unsigned char touch;
  float mxDst;
};








#endif /// GVETEX_H
