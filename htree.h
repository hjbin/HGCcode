//
// Created by huangjinbin on 2021/9/10.
//
#include <vector>
using namespace std;

#ifndef BOTTOM_UP_HTREE_H
#define BOTTOM_UP_HTREE_H
struct Node{
    int size;
    vector<Node *> child;
};

Node *newNode(int size)
{
    Node *temp=new Node;
    temp->size=size;
    return temp;
}
#endif //BOTTOM_UP_HTREE_H
