#ifndef GRAPH_H
#define GRAPH_H

#include <vector> 
#include <map>
#include <utility>

#include "Node.hpp"

class Graph {
    public:
        Graph(unsigned int k);
        void addSNV(std::string chrom, std::vector<unsigned int> pos, std::vector<std::string> ref, std::vector<std::string> alt);
        void addNode(std::string chrom, std::vector<unsigned int> pos, std::vector<std::string> bases, std::string suffix);
        void addEdge(std::string node1, std::string node2, int weight=1, std::string edgeColor="black", std::string direction="out");

        unsigned int getK();
        unsigned int getNodeCount();
        Node* getNode(std::string chrom, std::vector<unsigned int> pos, std::string suffix);
        Node* getNode(std::string nodeID);
        std::vector<uint> getOrderedSNVs();
        std::vector<Node*> getSNVNodes(std::string chrom, std::vector<unsigned int> pos);
        int getMaxWeightedEdge(const std::string& nodeID, const std::string& edgeColor, const std::string& direction);
        
        void populateGraph(std::string bcf_file_path);
        void connectNodeGroup(const std::vector<Node *>& nodes, std::string edgeColor);

        void printNodes();
        void printAdjList();
        void exportToDot(const std::string& filename);

        std::map<std::string, Node *> getNodes();
    private:
        static unsigned int nodeCount;
        unsigned int k;
        std::map<std::string, Node *> nodes;
        std::map<
            std::string, 
            std::pair<
                std::map< // stores out edges
                    std::string,
                    std::map<std::string, int>
                >, 
                std::map< // stores in edges
                    std::string,
                    std::map<std::string, int>
                >
            >
        > adjList;

        std::vector<std::string> orderedNodes;
        std::vector<uint> orderedSNVs;
};

#endif