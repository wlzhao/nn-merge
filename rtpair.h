#ifndef RTPAIR_H
#define RTPAIR_H


struct RTPair {
  float dst;
  unsigned int m, m2;
  unsigned int M;
  unsigned int id;
  unsigned char flag;

  static bool LtComp (const RTPair &a, const RTPair &b)
  {
    return (a.dst < b.dst);
  }

  static bool LgComp (const RTPair &a, const RTPair &b)
  {
    return (a.dst > b.dst);
  }

};


struct LRTPair
{
  public:
   float dst;
   unsigned long id;
   static bool LtComp(const LRTPair &a, const LRTPair &b)
   {
       return (a.dst > b.dst);
   }
};


#endif
