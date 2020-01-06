#ifndef NNITEM_H
#define NNITEM_H

struct NNItem {
    public:
    unsigned int Id;
    float dst;
    static int llComparer(const NNItem *pt1, const NNItem *pt2)
    {
         if(pt1->dst > pt2->dst)
         {
             return 0;
         }else if(pt1->dst < pt2->dst)
         {
             return 1;
         }else{
             return 1;
         }
     }
};


#endif // NNITEM_H
