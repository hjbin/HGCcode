//
// Created by huangjinbin on 2021/9/15.
//

#include <vector>
#include "hcluster.h"
using namespace std;

#ifndef BOTTOM_UP_UTILS_H
#define BOTTOM_UP_UTILS_H

struct Merge{
    int x,y;
};

struct Dendrogram{
    vector<Merge> e;
};

double modularity_by_level(vector<vector<int>> community)
{
    double mod=0.0;
    for (int i = 0; i < community.size(); ++i) {
        vector<int> cur_com=community[i];//process by community
        unordered_set<int> V;
        for (int k = 0; k < cur_com.size(); ++k) {//construct a vertex mapping
            V.insert(cur_com[k]);
        }

        int in_count=0;
        int out_count=0;
        for (int l = 0; l < cur_com.size(); ++l) {//bfs to count edge number
            int u=cur_com[l];
            out_count+=g[u].size();
            for (int k = 0; k < g[u].size(); ++k) {
                int v=g[u][k].u;
                if (V.find(v)!=V.end())//inner edge
                {
                    if (v>u)
                    {
                        in_count++;//count once
//                            out_count++;
                    }
                } else{//out edge
//                        out_count++;
                }
            }
        }
        double eii=(double)in_count/m*2;
        double ai=(double) out_count/m;
        mod+=eii-ai*ai;
    }
    return mod;
}

double singel_community_mod(vector<int> cur_com)
{
    double mod=0.0;
    unordered_set<int> V;
    for (int k = 0; k < cur_com.size(); ++k) {//construct a vertex mapping
        V.insert(cur_com[k]);
    }

    int in_count=0;
    int out_count=0;
    for (int l = 0; l < cur_com.size(); ++l) {//bfs to count edge number
        int u=cur_com[l];
        out_count+=g[u].size();
        if (cur_com.size()==1) break;// inner edges should be 0 if a community only contains 1 vertex
        for (int k = 0; k < g[u].size(); ++k) {
            int v=g[u][k].u;
            if (V.find(v)!=V.end())//inner edge
            {
                if (v>u)
                {
                    in_count++;//count once
//                            out_count++;
                }
            } else{//out edge
//                        out_count++;
            }
        }
    }
    double eii=(double)in_count/m*2;
    double ai=(double) out_count/m;
    mod=eii-ai*ai;


    return mod;
}


double modularity_by_dendrogram(Dendrogram d)
{
    vector<vector<int>> community;
    vector<double> com_mod;
    int layer=0;
    double total_mod;
    for (int j = 0; j < n; ++j) {
        if(nodemap[j])
        {
            vector<int> c;
            c.push_back(j);
            community.push_back(c);
            double mod_single=singel_community_mod(c);
            total_mod+=mod_single;
            com_mod.push_back(mod_single);

        }
    }
    layer++;

    int last_num=0,inc_layer=0;
    double total_inc=0.0;
    for (int i = 0; i < d.e.size(); ++i) {
        layer++;
        double layer_mod=0.0;
        Merge edge=d.e[i];
        int x=edge.x;
        int y=edge.y;

        //directly merge the x'th community with the y'th community
        vector<int> xcom=community[x];
        vector<int> ycom=community[y];

        //merge two community and add to the end of the list
        xcom.insert(xcom.end(),ycom.begin(),ycom.end());
        community.push_back(xcom);
        com_mod.push_back(singel_community_mod(xcom));//compute mod of new community

        //release the storage of the previous community
        community[x].clear();
        community[y].clear();

        int com_num=0;
        for (int j = 0; j < community.size(); ++j) {
            if (community[j].size()>0)
            {
                layer_mod+=com_mod[j];
                com_num++;
            }
        }
        if (i==1)
        {
            last_num=com_num;
        } else if (i>1 &&i<d.e.size()-1){
            total_inc+=(double)last_num-com_num;
//            cout<<"inc: "<<last_num-com_num<<endl;
            inc_layer++;
            last_num=com_num;
        }
        cout<<"modularity of level "<<layer<<": "<<layer_mod<<endl;
        total_mod+=layer_mod;

    }

    cout<<"Average Modularity: "<<total_mod/layer<<endl;
    cout<<"Average Community Quantity Increase: "<<total_inc/inc_layer<<endl;
    return total_mod/layer;
}


#endif //BOTTOM_UP_UTILS_H
