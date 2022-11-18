#include <iostream>
#include <vector>

using namespace std;

typedef struct Node{
    double LLR = 0;
    double INIT = 0;
    vector<int> neighbor;
    vector<double> OUT;
} node;

class Graph {
    public:
        node* VAR;
        node* CHK;
        Graph(int N, int M);
        ~Graph();
};

Graph::Graph(int N, int M){
    VAR = new node[N];
    CHK = new node[M];

    for(int i = 0; i < N; i++) VAR[i].OUT.assign(M, 0);
    for(int i = 0; i < M; i++) CHK[i].OUT.assign(N, 0);
}

Graph::~Graph(){
    delete VAR;
    delete CHK;
}