#ifndef GRAPH_H
#define GRAPH_H

#include<vector>
#include<string>
#include<iostream>

class Graph{
    public: 
        int n;
        std::vector<std::vector<int>> adj;
        std::vector<std::string> label;

        Graph(int n);
        void loadLabels(const std::string& filename);
        void addEdge(int u, int v);
        void print() const;
};
#endif