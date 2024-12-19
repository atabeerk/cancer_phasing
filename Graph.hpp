#ifndef GRAPH_H
#define GRAPH_H

#include <vector> 
#include <map>

#include "Node.hpp"

class Graph {
    public:
        Graph(unsigned int k);
        void addSNV(std::string chrom, std::vector<unsigned int> pos, std::vector<std::string> ref, std::vector<std::string> alt);
        void addNode(std::string chrom, std::vector<unsigned int> pos, std::vector<std::string> bases, std::string suffix);

        std::vector<Node*> getSNVNodes(std::string chrom, std::vector<unsigned int> pos);
        Node* getNode(std::string chrom, std::vector<unsigned int> pos, std::string suffix);
        Node* getNode(std::string nodeID);

        void addEdge(std::string node1, std::string node2, int weight=1);
        void connectNodeGroup(const std::vector<Node *>& nodes);

        unsigned int getNodeCount();
        void printNodes();
        void printAdjList();
        void exportToDot(const std::string& filename);

        unsigned int getK(); 
    private:
        static unsigned int nodeCount;
        unsigned int k;
        std::map<std::string, Node *> nodes;
        std::map<std::string, std::map<std::string, int>> adjList;
        std::vector<std::string> orderedNodes;
        std::vector<uint> orderedSNVs;
};

#endif