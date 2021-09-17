//
// Created by huangjinbin on 2021/9/14.
//

#ifndef BOTTOM_UP_DISJOINTSET_H
#define BOTTOM_UP_DISJOINTSET_H

struct SetTree{
    int father, rank, cost,count;
};
struct DisjointSet{
    SetTree *t;

    void init(int size)
    {
        t=new SetTree[size];
        for (int i = 0; i < size; ++i) {
            t[i].father=i;
            t[i].count=1;
            t[i].rank=0;
        }
    }

    int Find(int x)
    {
        while(t[x].father!=x)
        {
            x=t[x].father;
        }
        return x;
    }

    void Union(int x, int y)
    {
        int rx=Find(x);
        int ry=Find(y);

        if(rx==ry) return;
        if (t[rx].rank<t[ry].rank)
        {
            int temp=rx;
            rx=ry;
            ry=temp;
        }
        t[ry].father=rx;
        t[rx].count+=t[ry].count;
        if(t[rx].rank==t[ry].rank)
        {
            t[rx].rank++;
        }
    }

};
#endif //BOTTOM_UP_DISJOINTSET_H
