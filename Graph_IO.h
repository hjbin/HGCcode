#include "hcluster.h"
#include "Utils.h"
#include <iostream>
#include <sstream>
using namespace std;

//
// Created by huangjinbin on 2021/4/13.
//

#ifndef BOTTOM_UP_GRAPH_IO_H
#define BOTTOM_UP_GRAPH_IO_H

void Readin_Graph(char *str)
{
    int linecount=0;
    FILE *fin = fopen(str, "r");
    int x, y, z;
    if (skip)
    {
        cout<<"skip first line of input."<<endl;
        fscanf(fin, "%d %d", &x, &y);//skip the first line;
    }

    n = m = 0;
    dmax = 0;
    while(fscanf(fin, "%d %d", &x, &y)!=EOF)
    {
        n = max(n, x);
        n = max(n, y);
    }
    n++;
    nodemap=new bool [n];
    for (int k = 0; k < n; ++k) {
        nodemap[k]=false;
    }
    fclose(fin);
    g = new vector <endpoint> [n+1];

    //Vertex id starting from 1.
    fin = fopen(str, "r");
    if (skip)
    {

        fscanf(fin, "%d %d", &x, &y);//skip the first line;
    }
    endpoint a, b;
    Edgetype e;
    int eid = 0;

    while(fscanf(fin, "%d %d", &x, &y)!=EOF)
    {
        linecount++;
        if(linecount % 10000000==0)
        {
            cout<<"processed edge no: "<<linecount<<endl;
        }
        if (x==y){ continue;}
        if (y < x) swap(x, y);
        if (x < y)
        {
            nodemap[x]=true;
            nodemap[y]=true;
            if (eset.find(pair<int, int>(x, y)) == eset.end())
            {

                eset[pair<int, int>(x, y)] = eid;
                eset[pair<int, int>(y, x)] = eid;
                Edgetype ed;
                ed.x=x;
                ed.y=y;
                eg.push_back(ed);
                a.u = y;
                a.eid = eid;
                g[x].push_back(a);
                b.u = x;
                b.eid = eid;
                g[y].push_back(b);


                eid ++;
                m +=2;
            }
        }
        n = max(y, n);
    }

//    n++;

//    eg = new Edgetype [eid];
//    sup = new int [eid];
//    truss = new int [eid];

//    for (int i=0; i<n; i++)
//    {
//        int x= i;
//        for (int j=0; j<g[i].size(); j++)
//        {
//            int y = g[i][j].u;
//            int z = g[i][j].eid;
//            if (x< y)
//            {
//                e.x = x;
//                e.y = y;
//                eg.push_back(e);
////                eg[z] = e;
//                //printf("%d %d id=%d\n", x, y, z);
//            }
//
//        }
//    }


    for (int i=0; i<n; i++)
    {
        dmax = max(dmax, (int)g[i].size());
    }
    //	printf("n = %d, m = %d, dmax = %d, eid=%d\n", n, m, dmax, eid);
    fclose(fin);
}

void Readin_Paris(char *str)
{
    FILE *fin = fopen(str, "r");
    Communities.clear();
    vector<int> com;
    int x;
    int *interlabels=new int[n];
    int *lcs=new int[eg.size()];
    double dasgupta_cost=0.0;
    int ncount=0;
    for (int i = 0; i < eg.size(); ++i) {
        lcs[i]=eg.size();
    }
    for (int k = 0; k < n; ++k) {
        interlabels[k]=-1;
        if (nodemap[k])
        {
            ncount++;
        }
    }
    while (fscanf(fin, "%d", &x)!=EOF){
        if (x==-1)
        {
            Communities.push_back(com);
            com.clear();
        } else if (x==-2){
            for (int j = 0; j < Communities.size(); ++j) {//assign labels for nodes in one community layer
                vector<int> single_com=Communities[j];
                for (int k = 0; k < single_com.size(); ++k) {
                    interlabels[single_com[k]]=j;
                }
            }
            //Process by edge
            for (int l = 0; l < eg.size(); ++l) {
                Edgetype edg=eg[l];
                int x=edg.x;
                int y=edg.y;
                if (interlabels[x]==interlabels[y])
                {
                    if(lcs[l]>Communities[interlabels[x]].size())
                    {
                        lcs[l]=Communities[interlabels[x]].size();
                    }
                }
            }
//            cout<<Communities.size()<<endl;
//            HC.push_back(Communities);

            Communities.clear();
        } else{
            com.push_back(x);
        }
    }
    for (int k = 0; k < eg.size(); ++k) {
//        cout<<lcs[k]<<endl;
        dasgupta_cost+=(double)lcs[k];
    }

    printf("dasgupta_cost: %f\n", dasgupta_cost/(ncount*eg.size()));
}

void Readin_Paris_GT(char *str)
{
    FILE *fin = fopen(str, "r");
    Communities.clear();
    vector<int> com;
    vector<vector<int>> exact_com;
    int x;
    int *interlabels=new int[n];
    int *lcs=new int[eg.size()];
    bool visited_nodes=new bool[n];
    double dasgupta_cost=0.0;
    int ncount=0;
    for (int i = 0; i < eg.size(); ++i) {
        lcs[i]=eg.size();
    }
    for (int k = 0; k < n; ++k) {
        visited_nodes=false;
        interlabels[k]=-1;
        if (nodemap[k])
        {
            ncount++;
        }
    }
    while (fscanf(fin, "%d", &x)!=EOF){
        if (x==-1)
        {
            Communities.push_back(com);
            com.clear();
        } else if (x==-2){
            for (int j = 0; j < Communities.size(); ++j) {//assign labels for nodes in one community layer
                vector<int> single_com=Communities[j];
                if (single_com.size()>1)//only process merge
                {
                    for (int p=0;p<single_com.size();p++)
                    {

                    }
                }
            }
            //Process by edge
//            for (int l = 0; l < eg.size(); ++l) {
//                Edgetype edg=eg[l];
//                int x=edg.x;
//                int y=edg.y;
//                if (interlabels[x]==interlabels[y])
//                {
//                    if(lcs[l]>Communities[interlabels[x]].size())
//                    {
//                        lcs[l]=Communities[interlabels[x]].size();
//                    }
//                }
//            }
//            cout<<Communities.size()<<endl;
//            HC.push_back(Communities);

            Communities.clear();
        } else{
            com.push_back(x);
        }
    }
//    for (int k = 0; k < eg.size(); ++k) {
////        cout<<lcs[k]<<endl;
//        dasgupta_cost+=(double)lcs[k];
//    }
//
//    printf("dasgupta_cost: %f\n", dasgupta_cost/(ncount*eg.size()));
}


Dendrogram Input_Dendrogram(char *str)
{
    FILE *fin = fopen(str, "r");
    int x,y;
    Dendrogram d;
    while(fscanf(fin, "%d %d", &x, &y)!=EOF)
    {
        Merge me;
        me.x=x;
        me.y=y;
        d.e.push_back(me);
    }

    return d;
}

vector<vector<int>> Input_GTCom(char *str)
{
    string line;
    ifstream myfile;
    myfile.open (str);
    vector<vector<int>> gtComs;
    com_hash=new vector<int>[n];
    while (std::getline(myfile, line))
    {
        std::istringstream iss(line);
        int a;
        vector<int> gt_com;
        while (iss>>a)
        {
            gt_com.push_back(a);
            com_hash[a].push_back(gtComs.size());
        }
        sort(gt_com.begin(),gt_com.end());
        gtComs.push_back(gt_com);
//        if (gt_com.size()>3)
//        {
//            sort(gt_com.begin(),gt_com.end());
//            gtComs.push_back(gt_com);
//        }

    }
    return gtComs;

}
#endif //BOTTOM_UP_GRAPH_IO_H
