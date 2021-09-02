#include <iostream>
#include "Graph_IO.h"
#include <algorithm>
#include <math.h>
using namespace std;

void bottom_up()
{
    //    Readin_Graph("data/flights.txt");
    //SuperGraph sg=createSuperGraph();

    clock_t start_t, end_t;
    double otime=0.0;
    clock_t ostart,oend;

    start_t = clock();

    core_decompostion();
    SuperGraph sg;
    sg.init(g);
    sg.LabelProp();
    sg.retrieveCommunity();
    ostart=clock();
    //first layer
    vector<vector<int>> fl;
    for (int i = 0; i < n; ++i) {
        if(nodemap[i])
        {
            vector<int> c;
            c.push_back(i);
            fl.push_back(c);
        }
    }
    HC.push_back(fl);
    HC.push_back(Communities);
//    cout<<"Community Size: "<<Communities.size()<<endl;
//    cout<<"[";
//    for (int i = 0; i < Communities.size(); ++i) {
//        if (i==0){
//            cout<<"[";
//        } else{
//            cout<<",[";
//        }
//        for (int j = 0; j < Communities[i].size(); ++j) {
//            cout<<Communities[i][j];
//            if (j==Communities[i].size()-1){
//                cout<<"]";
//            } else{
//                cout<<",";
//            }
//        }
//    }
//    cout<<"]"<<endl;
    oend=clock();
    otime+=(double)(oend-ostart);
    sg.mergeSupernodes();
    while(sg.msize>0){
        sg.LabelProp();
        sg.retrieveCommunity();
        ostart=clock();
        HC.push_back(Communities);
//        cout<<"Community Size: "<<Communities.size()<<endl;
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
        oend=clock();
        otime+=(double)(oend-ostart);
        sg.mergeSupernodes();
    }



//    LabelPropagation(sg);
//    retrieveCommunity(sg);




    end_t = clock();
    double total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
    total_t=total_t-otime/CLOCKS_PER_SEC;
    printf("Total Time: =%.3lf\n", total_t);

//    outputHC();
}

void top_down_search()
{
    clock_t start_t, end_t;
    K=3;
    start_t=clock();
    core_decompostion();
    computeInformationGain();
    top_down();
    end_t=clock();
    double total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
    total_t=total_t-iotime/CLOCKS_PER_SEC;
    printf("Total Time: =%.3lf\n", total_t);

    //last layer
    vector<vector<int>> last_layer;
    for (int i = 0; i < n; ++i) {
        if (nodemap[i])
        {
            vector<int> c;
            c.push_back(i);
            last_layer.push_back(c);
        }
    }
    HC.push_back(last_layer);

    std::reverse(HC.begin(),HC.end());


//    outputHC();

}

int main(int argc, char* argv[]) {


    Readin_Graph(argv[1]);
    vector<vector<int>> gtCom=Input_GTCom(argv[2]);
//    for (int i = 0; i < 3; ++i) {
//        for (int j = 0; j < gtCom[i].size(); ++j) {
//            cout<<gtCom[i][j]<<" ";
//        }
//        cout<<endl;
//    }
//    Readin_Paris("data/.txt");
//    outputHC();
//    bottom_up();
    top_down_search();
    fast_dasgupta();
    hierarchical_fmeasure(gtCom);

//    core_decompostion();
//    std::cout << "Hello, World!" << std::endl;
    return 0;
}
