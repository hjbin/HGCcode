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

using namespace std;

int iter_T=10;

int Rand_Seed=10086;

int *labels;


int maxcore=0;

set<pair<double,int>> gain_sort;
vector<int> hierarchy;
vector<int> *core_bin;

vector<vector<int>> Communities;
vector<vector<vector<int>>> HC;

double iotime=0.0;


typedef struct endpoint
{
    int u;
    int eid;

}endpoint;

vector <endpoint> *g;
set<int> compset;
int n, m, dmax, maxsup, maxtruss;
int Tsize,K;


vector <int> *bin_sort;
vector<vector<int>> comm;
map<int,int> communitymapping;
int communityID;

bool* nodemap;

unordered_set<int> sparseEdges;

int test_num=0;
//flag[x]=true, the score of v has been computed.
bool *flag;
//conp[x] as the score of x
int *conp;
int *mv;
vector <endpoint> *inde;

int *coreness;
vector<int> qnodes;
int *nodetruss;

//the sum of Top-k values
int total_comp =0;

set<int> removedEdges;

typedef struct Edgetype
{
    int x, y,truss;
}Edgetype;


vector<Edgetype> eg;
int *sup, *truss;

map<int,vector<endpoint>> *t;

map <pair<int, int>, int > eset;
vector <int> *L;

string graphname;

bool *del;

typedef struct Score
{
    int id;
    int score;
}Score;

bool operator< (const Score& t1, const Score& t2) {
    return t1.score > t2.score;
}

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

typedef struct AdvancedInde
{
    //vector<Edgetype> superedges;
    map<int,SuperNode> supernodes;
    map<int,vector<int>> nodebin;
    map<int,vector<Edgetype>> edgebin;
    int knmax;
    int kmmax;
    //vector<int> internal_nodes;
}AdvancedInde;

AdvancedInde *adi;


typedef struct Index
{
    double upper_bound;
    map <int, Setree> t;

    int change;

    void init(int id) {
        t[0].father = 0;
        t[0].rank = -1;
        t[0].count = 0;
        t[0].cost = 0;

        for (int i=0; i<g[id].size(); i++)
        {
            int x = g[id][i].u;
            t[x].father = x;
            t[x].rank = 0;
            t[x].count = 1;
            t[x].cost = min(g[x].size(), g[id].size());
        }

    }

    void initSubgraph(int nid){
        t[0].father = 0;
        t[0].rank = -1;
        t[0].count = 0;
        t[0].cost = 0;
        for (int i = 1; i < nid; ++i) {
            t[i].father = i;
            t[i].rank = 0;
            t[i].count = 1;
            t[i].cost = 0;
        }
    }

    int Find_Father(int x)
    {
        if (x != t[x].father)
        {
            t[x].father = Find_Father(t[x].father);
            return t[x].father;
        }
        return t[x].father;
    }
    void Union(int x,int y)
    {
        int fx, fy;
        fx = Find_Father(x);
        fy = Find_Father(y);

        if (fx == fy) return ;
        if (t[fx].rank > t[fy].rank)
        {
            t[fy].father = fx;
            t[fx].count += t[fy].count;
            t[fx].cost += t[fy].cost;
        }
        else
        {
            t[fx].father = fy;
            t[fy].count += t[fx].count;
            t[fy].cost += t[fx].cost;

            if (t[fx].rank == t[fy].rank)
            {
                t[fy].rank++;
            }
        }
        change = 1;
    }
    void Single_Father(int x)
    {
        t[x].father = 0;
        t[0].count++;
        change = 1;
    }
}Index;

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
                      labels[j]=now_label[j]=get_maxnum_in_vec(neighbor_labels);
//                    labels[j]=now_label[j]=get_maxnum_in_vec_with_core(neighbor_labels,neighbor_core);
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


//Union_Find_Isolate data structure(AB)
Index *Set;






void core_decompostion()
{
    //cout<<g[vertex].size()<<endl;

    clock_t start_t, end_t;
    maxsup = 0;
    maxtruss = 0;
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
                    //neighbor_labels.push_back(pre_label[sgraph.g[j][k].u]);//sync
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

void hierarchical_fmeasure(vector<vector<int>> gtCom)
{
    cout<<"gt com size: "<<gtCom.size()<<endl;
    int hsum=(1+HC.size())*HC.size()/2;
    double total_fk=0.0;
    for (int i = 1; i < HC.size(); ++i) {

        int layer=HC.size()-i;
        cout<<"fscore layer "<<layer<<endl;
        double alpha=(double)layer/hsum;
        double layer_fk=0.0;
        for (int j = 0; j < HC[i].size(); ++j) {
            vector<int> C=HC[i][j];
            sort(C.begin(),C.end());
            layer_fk+=alpha*(F_k(1,C,gtCom)/HC[i].size());
        }
        total_fk+=layer_fk/(HC.size()-1);
    }
    cout<<"Hierarchical Fmeasure: "<<total_fk<<endl;
}





//int select_max_labels(vector<int> v_label){
//    int label=0;
//
//
//    return label;
//}

#endif //BOTTOM_UP_HCLUSTER_H
