//
// Created by huangjinbin on 2021/4/13.
//

#ifndef BOTTOM_UP_HCLUSTER_H
#define BOTTOM_UP_HCLUSTER_H

#include <fstream>
#include <unordered_set>
#include <vector>
#include <queue>
#include <set>
#include <map>
#include <unordered_map>
#include <algorithm>
#include "htree.h"

using namespace std;

int iter_T=10;

int Rand_Seed=10086;

int *labels;
bool skip= false;
bool hlp=false;

int maxcore=0;

set<pair<double,int>> gain_sort;
vector<int> hierarchy;
vector<int> *core_bin;

vector<vector<int>> Communities;
vector<vector<vector<int>>> HC;
vector<int> *com_hash;

double iotime=0.0;


typedef struct endpoint
{
    int u;
    int eid;

}endpoint;

vector <endpoint> *g;
set<int> compset;
int n, m, dmax;
int Tsize,K;


vector <int> *bin_sort;
vector<vector<int>> comm;
map<int,int> communitymapping;
int communityID;

bool* nodemap;

unordered_set<int> sparseEdges;


vector <endpoint> *inde;

int *coreness;


Node *root;

set<int> removedEdges;

typedef struct Edgetype
{
    int x, y;
}Edgetype;


vector<Edgetype> eg;

map<int,vector<endpoint>> *t;

map <pair<int, int>, int > eset;
vector <int> *L;

string graphname;

bool *del;

typedef struct Setree
{
    int father, rank, count;
    int cost;
}Setree;

typedef struct SuperNode
{
    int root_id;
    vector<int> internal_nodes;
    vector<int> nb;
    int coreness;
}SuperNode;

int Random(int l, int r){srand(Rand_Seed); return rand() % (r - l + 1) + l;}

struct SuperGraph{
    vector<endpoint> *g;
    vector<Edgetype> edge;
    int size;
    int msize;
    int *core;
    vector<vector<int>> scommunity;
    vector<SuperNode> nlist;


    void init(vector<endpoint> *gt){
        for (int i = 0; i < n; ++i) {
            SuperNode sn;
            sn.internal_nodes.push_back(i);
            sn.root_id=i;
            sn.coreness=coreness[i];
            for (int j = 0; j < gt[i].size(); ++j) {
                sn.nb.push_back(gt[i][j].u);
            }
            nlist.push_back(sn);
        }
        core=coreness;
        size=n;
        msize=m;
        edge=eg;
    }

    void LabelProp(){
        labels=new int[size];
        vector<int> pre_label,now_label;
        for (int i = 0; i < size; ++i) {
            labels[i]=i;
            now_label.push_back(labels[i]);
        }
        for (int i = 0; i < iter_T; ++i) {

            pre_label.assign(now_label.begin(), now_label.end());
            for (int j = 0; j < size; ++j) {
                vector<int> neighbor_labels;
                vector<int> neighbor_core;
                for (int k = 0; k < nlist[j].nb.size(); ++k) {
//                    neighbor_labels.push_back(pre_label[nlist[j].nb[k]]);//sync
                    neighbor_labels.push_back(now_label[nlist[j].nb[k]]);//async
                    neighbor_core.push_back(core[nlist[j].nb[k]]);
                }
                if (nlist[j].nb.size()>0){
                    if(hlp){
                        labels[j]=now_label[j]=get_maxnum_in_vec(neighbor_labels);
                    }else{
                    labels[j]=now_label[j]=get_maxnum_in_vec_with_core(neighbor_labels,neighbor_core);
                    }
                }
            }
            if (pre_label==now_label) break;
        }
    }

    void mergeSupernodes(){
        vector<SuperNode> nodes;
        vector<Edgetype> tmp_edge;
        map<pair<int,int>,int> emap;
        map<int,int> nodemap;
        int mid=0;
        int id=size+1;
        int nid=0;
        for (int i = 0; i < scommunity.size(); ++i) {
            SuperNode newNode;
            newNode.root_id=nid;
            int maxcore=0;
            for (int j = 0; j < scommunity[i].size(); ++j) {
                int rid=scommunity[i][j];
                if(maxcore<nlist[rid].coreness){
                    maxcore=nlist[rid].coreness;
                }
                nlist[rid].root_id=id;
                newNode.internal_nodes.insert(newNode.internal_nodes.end(),nlist[rid].internal_nodes.begin(),nlist[rid].internal_nodes.end());
                nodemap[id]=nid;
            }
            newNode.coreness=maxcore;
            nodes.push_back(newNode);
            id++;
            nid++;
        }
        for (int k = 0; k < edge.size(); ++k) {
            Edgetype ed=edge[k];
            int x=ed.x;
            int y=ed.y;
            if (nlist[x].root_id!=nlist[y].root_id){
                int x1=nodemap[nlist[x].root_id];
                int y1=nodemap[nlist[y].root_id];
                if (emap.find(pair<int,int>(x1,y1))==emap.end()){
                    emap[pair<int,int>(x1,y1)]=mid;
                    emap[pair<int,int>(y1,x1)]=mid;
                    nodes[x1].nb.push_back(y1);
                    nodes[y1].nb.push_back(x1);
                    Edgetype superedge;
                    superedge.x=x1;
                    superedge.y=y1;
                    tmp_edge.push_back(superedge);
                    mid++;
                }

            }
        }

        nlist.clear();
        nlist=nodes;
        edge.clear();
        edge=tmp_edge;
        size=nlist.size();
        msize=edge.size();
        emap.clear();
        cout<<"super graph size (n/m): "<<"("<<size<<"/"<<msize<<")"<<endl;

    }

    void retrieveCommunity()
    {
//    int label_size = sizeof(labels)/sizeof(labels[0]);
        vector<int> *comm=new vector<int>[size];
        vector<int> *scomm=new vector<int>[size];
        scommunity.clear();
        Communities.clear();
        for (int i = 0; i < size; ++i) {
            bool nullflag=false;
//        cout<<"node id: "<<i<<", label: "<<labels[i]<<endl;
            for (int j = 0; j < nlist[i].internal_nodes.size(); ++j) {
                int real_i=nlist[i].internal_nodes[j];
                if (!nodemap[real_i]){
                    continue;
                }
                comm[labels[i]].push_back(real_i);
                nullflag=true;
            }
            if (nullflag){
                scomm[labels[i]].push_back(i);
            }

        }
        for (int j = 0; j < size; ++j) {
            if (comm[j].size()>0){
                Communities.push_back(comm[j]);
                scommunity.push_back(scomm[j]);
            }
        }
        delete [] labels;
        delete [] comm;
        delete [] scomm;
    }

    int get_maxnum_in_vec_with_core(vector <int> vec,vector<int> core){
        map <int, int> mp; // need to optimize
        int Max = 0;
        vector <int> max_vec;
        for (int i = 0; i < vec.size(); i++){
            mp[vec[i]]+=core[i];
            if (mp[vec[i]] > Max){
                Max = mp[vec[i]];
                max_vec.clear();
            }
            if (mp[vec[i]] == Max) max_vec.push_back(vec[i]);
        }
        //cout << max_vec.size() << endl;
        return max_vec[Random(0, max_vec.size() - 1)];
    }

    int get_maxnum_in_vec(vector <int> vec){
        map <int, int> mp; // need to optimize
        int Max = 0;
        vector <int> max_vec;
        for (int i = 0; i < vec.size(); i++){
            mp[vec[i]]++;
            if (mp[vec[i]] > Max){
                Max = mp[vec[i]];
                max_vec.clear();
            }
            if (mp[vec[i]] == Max) max_vec.push_back(vec[i]);
        }
        //cout << max_vec.size() << endl;
        return max_vec[Random(0, max_vec.size() - 1)];
    }

    void mergeNodes(){
        vector<SuperNode> newNode;
        unordered_map<int,SuperNode> mp;
        vector<Edgetype> newEdge;
        int sn_id=0;
        for (int i = 0; i < msize; ++i) {
            Edgetype ed=edge[i];
            int x=ed.x;
            int y=ed.y;

            if (labels[x]==labels[y]){
                SuperNode sx=mp[x];
                SuperNode sy=mp[y];

                if (sx.root_id!=sy.root_id){

                }

            } else{
                newEdge.push_back(ed);
            }



        }
    }
};








void core_decompostion()
{
    //cout<<g[vertex].size()<<endl;

    clock_t start_t, end_t;
    L = new vector <int> [n+1];



//    sup=new int[m];
//    truss=new int[m];

    int *degg=new int [n+1];
    int maxdeg=0;

//    for (int i1 = 0; i1 < edgeid; ++i1) {
//        sup[i1]=0;
//    }

    removedEdges.clear();


    start_t = clock();

    for (int i=0; i<n+1; i++)
    {
        int x = i;
        degg[x]=g[x].size();
        if (maxdeg<degg[x])
        {
            maxdeg=degg[x];
        }
        L[degg[x]].push_back(x);

    }

//    for (int i=0; i<g[vertex].size(); i++)
//    {
//        int x = g[vertex][i].u;
//        for (int j=0; j<g[x].size(); j++)
//        {
//            int y = g[x][j].u;
//            if (eset.find(pair<int, int>(vertex, y)) == eset.end())
//            {
//                continue;
//                //cout<<u<<"  "<<v<<"  "<<w<<" "<<sup[eid]<<endl;
//            }
//            if (x <y )
//            {
//                int eid = g[x][j].eid;
//                sup[eid] = 0;
//                int u, v;
//                u = y;
//                v = x;
//                if (g[x].size() <g[y].size())
//                {
//                    u  = x;
//                    v = y;
//                }
//                for (int k=0; k<g[u].size(); k++)
//                {
//                    int w = g[u][k].u;
//                    if (eset.find(pair<int, int>(v, w)) != eset.end())
//                    {
//
//                        if (eset.find(pair<int, int>(vertex, u)) != eset.end() && eset.find(pair<int, int>(vertex, v)) != eset.end() && eset.find(pair<int, int>(vertex, w)) != eset.end() )
//                        {
//                            sup[eid]++;
//                            //cout<<u<<"  "<<v<<"  "<<w<<" "<<sup[eid]<<endl;
//                        }
//                        //cout<<eid<<endl;
//                    }
//                }
//                //printf("%d %d id=%d, sup=%d\n", x, y, eid, sup[eid]);
//                if ((sup[eid]) > maxsup)
//                {
//                    maxsup = sup[eid];
//                }
//                L[sup[eid]].push_back(eid);
//            }
//        }
//    }
    end_t = clock();


    double total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
    //printf("Triangle Count: =%.3lf\n", total_t);

//    for (int i=0;i<m/2;i++)
//    {
//        cout<<sup[i]<<endl;
//    }

    start_t = clock();
    del = new bool [n+1];
    coreness=new int[n+1];
    for (int i=0; i<n+1; i++)
    {
        del[i] = false;
    }
    int kmax=0;
//	vector<Edgetype> *gx=new vector<Edgetype> [maxsup+3];
    int md=maxdeg;
//	if(maxdeg>Tsize){
//		maxdeg=Tsize;
//	}
    for (int k=0; k<=maxdeg; k++)
    {
        for (int i=0; i<L[k].size(); i++)
        {
            int x=L[k][i];
            if (del[x]==true)
            {
                continue;
            }
            coreness[x]=k;
            if (maxcore<coreness[x])
            {
                maxcore=coreness[x];
            }
            del[x]=true;
            for (int j = 0; j < g[x].size(); ++j) {
                int y=g[x][j].u;

                if (del[y]==true)
                {
                    continue;
                }

                int realx=x;
                int realy=y;

                int delEid=eset[pair<int,int>(realx,realy)];

                if (coreness[x]<Tsize)
                {
                    removedEdges.insert(delEid);
                }

                int temp=y;

                if (degg[temp]>k)
                {
                    degg[temp]--;
                    L[degg[temp]].push_back(temp);
                }


            }
            //printf("%d %d id=%d, truss=%d\n", x, y, eid, truss[eid]);
        }
    }


    end_t = clock();
    total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
    printf("Core Decomposition: =%.3lf\n", total_t);
    //cout<<"id="<<vertex<<", ";
    printf("maxdegree=%d, maxcore=%d \n", md, maxcore);
    delete [] L;
    delete [] inde;
    delete [] del;
//    delete [] sup;
//    delete [] truss;
    delete [] degg;
}


int get_maxnum_in_vec_with_core(vector <int> vec,vector<int> core){
    map <int, int> mp; // need to optimize
    int Max = 0;
    vector <int> max_vec;
    for (int i = 0; i < vec.size(); i++){
        mp[vec[i]]+=core[i];
        if (mp[vec[i]] > Max){
            Max = mp[vec[i]];
            max_vec.clear();
        }
        if (mp[vec[i]] == Max) max_vec.push_back(vec[i]);
    }
    //cout << max_vec.size() << endl;
    return max_vec[Random(0, max_vec.size() - 1)];
}

int get_maxnum_in_vec(vector <int> vec){
    map <int, int> mp; // need to optimize
    int Max = 0;
    vector <int> max_vec;
    for (int i = 0; i < vec.size(); i++){
        mp[vec[i]]++;
        if (mp[vec[i]] > Max){
            Max = mp[vec[i]];
            max_vec.clear();
        }
        if (mp[vec[i]] == Max) max_vec.push_back(vec[i]);
    }
    //cout << max_vec.size() << endl;
    return max_vec[Random(0, max_vec.size() - 1)];
}

void LabelPropagation(SuperGraph sgraph){
    labels=new int[sgraph.size];
    vector<int> pre_label,now_label;
    for (int i = 0; i < sgraph.size; ++i) {
        labels[i]=i;
        now_label.push_back(labels[i]);
    }
    for (int i = 0; i < iter_T; ++i) {

        pre_label.assign(now_label.begin(), now_label.end());
        for (int j = 0; j < sgraph.size; ++j) {
            vector<int> neighbor_labels;
            vector<int> neighbor_core;
            for (int k = 0; k < sgraph.g[j].size(); ++k) {
                //neighbor_labels.push_back(pre_label[sgraph.g[j][k].u]);//sync
                neighbor_labels.push_back(now_label[sgraph.g[j][k].u]);//async
                neighbor_core.push_back(coreness[sgraph.g[j][k].u]);
            }
            if (sgraph.g[j].size()>0){
                labels[j]=now_label[j]=get_maxnum_in_vec(neighbor_labels);
//                labels[j]=now_label[j]=get_maxnum_in_vec_with_core(neighbor_labels,neighbor_core);
            }
        }
        if (pre_label==now_label) break;
    }

}

void computeInformationGain(){
    int coresize=0;
    core_bin=new vector<int>[maxcore+1];
    for (int i = 0; i < n; ++i) {
        core_bin[coreness[i]].push_back(i);
    }

    for (int j = maxcore; j >=1; j--)
    {
        if (j==maxcore)
        {
            double info_gain=0.0;
            coresize+=core_bin[maxcore].size();
            gain_sort.insert(pair<double,int>(info_gain,j));
        } else{
            coresize+=core_bin[j].size();
            double info_gain=(double)core_bin[j].size()/coresize;
            gain_sort.insert(pair<double,int>(info_gain,j));
//            cout<<"bin size: "<<core_bin[j].size()<<", core size: "<<coresize<<endl;
//            cout<<"gain: "<<info_gain<<endl;
        }
    }

//    for (int l = 0; l < maxcore+1; ++l) {
//        double info_gain=double(H[l].core.size()-H[l+1].core.size())/H[l].core.size();
////        cout<<info_gain<<","<<l<<endl;
//        gain_sort.insert(pair<double,int>(info_gain,l));
//    }


    int count=0;
    if (K>maxcore)
    {
        K=maxcore;
    }
    for (auto it=gain_sort.rbegin();it!=gain_sort.rend();it++) {
        if (count<K){
//            cout<<it->first<<endl;
            hierarchy.push_back(it->second);
        } else{
            break;
        }
        count++;
    }

}

void retrieveCommunity(int *labels)
{
    Communities.clear();
//    int label_size = sizeof(labels)/sizeof(labels[0]);
    vector<int> *comm=new vector<int>[n];
    for (int i = 0; i < n; ++i) {
//        cout<<"node id: "<<i<<", label: "<<labels[i]<<endl;
        if (!nodemap[i]){
            continue;
        }
        comm[labels[i]].push_back(i);
    }
    for (int j = 0; j < n; ++j) {
        if (comm[j].size()>0){
            Communities.push_back(comm[j]);
        }
    }
    if (HC.size()>0)
    {
        if (HC[HC.size()-1].size()!=Communities.size())
        {
            HC.push_back(Communities);
        }
    } else{
        HC.push_back(Communities);
    }
}

void output()
{
    clock_t ostart,oend;
    ostart=clock();
    cout<<"Community Size: "<<Communities.size()<<endl;
    cout<<"[";
    for (int i = 0; i < Communities.size(); ++i) {
        if (i==0){
            cout<<"[";
        } else{
            cout<<",[";
        }
        for (int j = 0; j < Communities[i].size(); ++j) {
            cout<<Communities[i][j];
            if (j==Communities[i].size()-1){
                cout<<"]";
            } else{
                cout<<",";
            }
        }
    }
    cout<<"]"<<endl;
    oend=clock();
    iotime+=(double)(oend-ostart);
}

void outputHC()
{
    for (int x = 0; x < HC.size(); ++x) {
        vector<vector<int>> Communities=HC[x];
        cout<<"Community Size: "<<Communities.size()<<endl;
//        for (int k = 0; k < Communities.size(); ++k) {
//            cout<<Communities[k].size()<<",";
//        }
//        cout<<endl;
//        cout<<"[";
//        for (int i = 0; i < Communities.size(); ++i) {
//            if (i==0){
//                cout<<"[";
//            } else{
//                cout<<",[";
//            }
//            for (int j = 0; j < Communities[i].size(); ++j) {
//                cout<<Communities[i][j];
//                if (j==Communities[i].size()-1){
//                    cout<<"]";
//                } else{
//                    cout<<",";
//                }
//            }
//        }
//        cout<<"]"<<endl;
    }
}



void top_down(){
    int *last_label=new int[n];
    int *cur_label=new int[n];
    vector<int> onlyCom;
    for (int i = 0; i < n; ++i) {
        if(nodemap[i])
        {
            onlyCom.push_back(i);
        }
        last_label[i]=1;
//        cur_label[i]=-1;
    }
    Communities.push_back(onlyCom);
    HC.push_back(Communities);
    sort(hierarchy.begin(),hierarchy.end());
    for (int i = hierarchy.size()-1; i >=0; --i) {
        int core_num=hierarchy[i];
        cout<<"layer k: "<<core_num<<endl;
        //select seeds from k-core
        vector<int> seeds;
        for (int j = core_num; j <= maxcore; ++j) {
            if (core_bin[j].size()==0){
                continue;
            }
            for (int k = 0; k < core_bin[j].size(); ++k) {
                seeds.push_back(core_bin[j][k]);
            }
        }

        //init label
        vector<int> now_label,pre_label;
        for (int p = 0; p < n; ++p) {
            cur_label[p]=-1;
        }
        for (int l = 0; l < seeds.size(); ++l)
        {
            if (i==0){//generate new label;
                cur_label[seeds[l]]=seeds[l];//ensure no duplicated label in lower hierarchies
                now_label.push_back(seeds[l]);
            } else{
                if(coreness[seeds[l]]>core_num)
                {
                    cur_label[seeds[l]]=last_label[seeds[l]];
                } else{
                    cur_label[seeds[l]]=seeds[l];
                }
                now_label.push_back(cur_label[seeds[l]]);

            }

        }

//        for (int k = 0; k < n; ++k) {
////            last_label[k]=cur_label[k];
//            cout<<cur_label[k]<<",";
//        }
//        cout<<endl;

        //perform LPA
        for (int p = 0; p < iter_T; ++p)
        {
            pre_label.assign(now_label.begin(),now_label.end());
            for (int j = 0; j < seeds.size(); ++j)
            {
                int x=seeds[j];
                vector<int> neighbor_labels;
                for (int k = 0; k < g[x].size(); ++k) {
                    int y=g[x][k].u;

                    if (last_label[x]!=last_label[y] || coreness[y]<core_num)//connectivity check
                    {
                        continue;
                    }
//                    neighbor_labels.push_back(pre_label[j]);//sync
                    neighbor_labels.push_back(cur_label[y]);//async
                }
                if (neighbor_labels.size()>0){
                    cur_label[x]=now_label[j]=get_maxnum_in_vec(neighbor_labels);
//                labels[j]=now_label[j]=get_maxnum_in_vec_with_core(neighbor_labels,neighbor_core);
                }
            }
            if (pre_label==now_label) break;
        }

//        for (int k = 0; k < n; ++k) {
////            last_label[k]=cur_label[k];
//            cout<<cur_label[k]<<",";
//        }
//        cout<<endl;

        //bfs assign label
        queue<int> q;
        for (int j = 0; j < n; ++j) {
            if(cur_label[j]==-1) continue;
            q.push(j);
        }

        while (!q.empty())
        {
            int x=q.front();
            q.pop();
            for (int j = 0; j < g[x].size(); ++j) {
                int y=g[x][j].u;
                if (cur_label[y]!=-1 || last_label[x]!=last_label[y])//connectivity check
                {
                    continue;
                }
                vector<int> neigbor_label;
//                neigbor_label.push_back(x);
                for (int k = 0; k < g[y].size(); ++k) {
                    int z=g[y][k].u;
                    if (cur_label[z]==-1 || last_label[y]!=last_label[z])//ensure neighbor has label and connected in last layer
                    {
                        continue;
                    }
                    neigbor_label.push_back(cur_label[z]);
                }
                cur_label[y]=get_maxnum_in_vec(neigbor_label);
                q.push(y);
            }
        }
        for (int k = 0; k < n; ++k) {
            last_label[k]=cur_label[k];
//            cout<<cur_label[k]<<",";
        }
//        cout<<endl;
        retrieveCommunity(cur_label);
//        output();
    }

    delete[] cur_label;
    delete[] last_label;
}







SuperGraph createSuperGraph()
{
    SuperGraph superGraph;
    superGraph.g=g;
    superGraph.size=n;
    return superGraph;
}


//complexity O(mh)
void fast_dasgupta()
{
    int *interlabels=new int[n];
    int *lcs=new int[eg.size()];
    double dasgupta_cost=0.0;
    int ncount=0;
    for (int k = 0; k < n; ++k) {
        interlabels[k]=-1;
        if (nodemap[k])
        {
            ncount++;
        }
    }
    for (int i = HC.size()-1; i >=0 ; --i) {
        vector<vector<int>> comm=HC[i];
        for (int j = 0; j < comm.size(); ++j) {//assign labels for nodes in one community layer
            vector<int> single_com=comm[j];
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
                lcs[l]=comm[interlabels[x]].size();
            }
        }
    }
    for (int k = 0; k < eg.size(); ++k) {
        dasgupta_cost+=(double)lcs[k];
    }

    printf("dasgupta_cost: %f\n", dasgupta_cost/m/ncount);

}

void lca_cost(int L=1)
{
    int *interlabels=new int[n];
    int *lcs=new int[eg.size()];
    double lca_cost=0.0;
    int ncount=0;
    for (int k = 0; k < n; ++k) {
        interlabels[k]=-1;
        if (nodemap[k])
        {
            ncount++;
        }
    }
    for (int i = HC.size()-1; i >=0 ; --i) {
        vector<vector<int>> comm=HC[i];
        for (int j = 0; j < comm.size(); ++j) {//assign labels for nodes in one community layer
            vector<int> single_com=comm[j];
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
                lcs[l]=i;
            }
        }
    }

    for (int k = 0; k < eg.size(); ++k) {
        lca_cost+=(double)lcs[k];
    }
    //cout<<"edge:"<<m<<" "<<eg.size()<<"\n";
    printf("lca_cost: %f\n", lca_cost/eg.size()/HC.size());

}




void treeSamplingDivergence(vector<vector<vector<int>>> L){

    double tsd=0.0;
    for (int l = 0; l < eg.size(); ++l) {
        int v=eg[l].x;
        int u=eg[l].y;
        int csize=n;
        for (int i1 = L.size(); i1 >= 0; i1--){
            int i=hierarchy[i1];
            bool vfind= false;
            bool ufind= false;

//            for (int j = 1; j <= L[i].size(); j++){
////            printf("Label %d: ", j);
//                for (int k = 0; k < L[i].community[j].size(); k++){
//                    if(L[i].community[j][k]==v){
//                        vfind= true;
//                    }
//                    if(L[i].community[j][k]==u){
//                        ufind= true;
//                    }
//                }
//                if (ufind && vfind){
//                    if (L[i].community[j].size()<csize){
//                        csize=L[i].community[j].size();
//                    }
//                    break;
//                }
//
//            }
        }
        tsd+=((double)1/m)*((double)csize/n);


    }

    printf("dasgupta_cost: %f\n", tsd);
}

vector<int> set_intersect(vector<int> set1, vector<int> set2)
{
    vector<int> intersection;

    std::set_intersection(set1.begin(),set1.end(),set2.begin(),set2.end(),std::back_inserter(intersection));

    return intersection;
}

double F_k(int k,vector<int> C1,vector<vector<int>> C2)
{
    double max_fk=0.0;
    for (int i = 0; i < C2.size(); ++i) {
        vector<int> cstar=C2[i];
        vector<int> inter=set_intersect(C1,cstar);
        double precision=(double)inter.size()/C2.size();
        double recall=(double)inter.size()/cstar.size();
        double Fk=(k+1)*precision*recall/(k*precision+recall);
        if (max_fk<Fk)
        {
            max_fk=Fk;
        }
    }
    return max_fk;
}

double F_k_fast(int k,vector<int> C1,vector<vector<int>> gt)
{
    double max_fk=0.0;
    unordered_set<int> s;
    vector<vector<int >> C2;
    for (int i = 0; i < C1.size(); ++i) {
        int x=C1[i];
        for (int j = 0; j < com_hash[x].size(); ++j) {
            if (s.find(com_hash[x][j])!=s.end())
            {
                continue;
            }
            s.insert(com_hash[x][j]);
            C2.push_back(gt[com_hash[x][j]]);
        }
    }


    for (int i = 0; i < C2.size(); ++i) {
        vector<int> cstar=C2[i];
        vector<int> inter=set_intersect(C1,cstar);
        double precision=(double)inter.size()/C1.size();
        double recall=(double)inter.size()/cstar.size();
        double Fk=(k+1)*precision*recall/(k*precision+recall);
        if (max_fk<Fk)
        {
            max_fk=Fk;
        }
    }
    return max_fk;
}

void hierarchical_fmeasure(vector<vector<int>> gtCom)
{
    cout<<"gt com size: "<<gtCom.size()<<endl;
    int hsum=(1+HC.size())*HC.size()/2;
    double total_fk=0.0;
    for (int i = 0; i < HC.size(); ++i) {

        int layer=HC.size()-i;
        cout<<"fscore layer "<<layer<<endl;
        double alpha=(double)layer/hsum;
        double layer_fk=0.0;
        for (int j = 0; j < HC[i].size(); ++j) {
            vector<int> C=HC[i][j];
            sort(C.begin(),C.end());
            layer_fk+=alpha*(F_k_fast(1,C,gtCom)/HC[i].size());
        }
        total_fk+=layer_fk/(HC.size());
    }
    cout<<"Hierarchical Fmeasure: "<<total_fk<<endl;
}

void hierarchical_gain()
{
    double total=0.0;
    int layer=0;
    for (int i = 1; i < HC.size()-2; ++i) {
        double increase=(double)HC[i].size()-HC[i+1].size();
        layer++;
//        cout<<increase<<endl;
        total+=increase;
    }
    cout<<"Average Community Quantity Increase: "<<total/layer<<endl;
}

void average_modularity()
{
    double total_mod=0.0;
    for (int i = 0; i < HC.size(); ++i) {
        vector<vector<int>> H=HC[i];
        double mod=0.0;

        for (int j = 0; j < H.size(); ++j) {
            vector<int> cur_com=H[j];//process by community
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
        cout<<"modularity of level "<<i<<": "<<mod<<endl;
        total_mod+=mod;
    }

    cout<<"average modularity: "<<total_mod/HC.size()<<endl;
}



//int select_max_labels(vector<int> v_label){
//    int label=0;
//
//
//    return label;
//}

#endif //BOTTOM_UP_HCLUSTER_H
