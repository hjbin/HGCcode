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

double communitySizeIncrease(Node *root)
{
    queue<Node *> q;
    queue<Node *> nextQ;
    q.push(root);
    int level=0;
    int sizeInc=0;
    while(!q.empty())
    {
        Node *cur=q.front();
        q.pop();
        int maxSize=0;
        for (int i = 0; i < cur->child.size(); ++i) {
            nextQ.push(cur->child[i]);
            if (maxSize<cur->child[i]->size)
            {
                maxSize=cur->child[i]->size;
            }
        }
        sizeInc+=cur->size-maxSize;
        if (q.empty())
        {
            q=nextQ;
            level++;
            nextQ=queue<Node *>();
        }
    }

    double avgInc=(double)sizeInc/level;
    return avgInc;
}

#endif //BOTTOM_UP_HTREE_H
