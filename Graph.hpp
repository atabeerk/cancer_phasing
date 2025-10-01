#ifndef GRAPH_H
#define GRAPH_H

#include <vector> 
#include <map>
#include <utility>

#include "Node.hpp"
#include "config.hpp"

class Graph {
    public:
        Graph(unsigned int k);
        void addSNV(const std::string& chrom, const std::vector<unsigned int>& pos, const std::vector<std::string>& ref, const std::vector<std::string>& alt);
        void addNode(std::vector<std::string>& layer, std::string chrom, std::vector<unsigned int> pos, std::vector<std::string> bases, std::string suffix);
        void addEdge(std::string node1, std::string node2, int weight=1, std::string edgeColor="black", std::string direction="out");

        unsigned int getK();
        unsigned int getNodeCount();
        Node* getNode(std::string chrom, std::vector<unsigned int> pos, std::string suffix);
        Node* getNode(std::string nodeID);
        std::vector<uint> getOrderedSNVs();
        std::vector<Node*> getSNVNodes(const std::string& chrom, const std::vector<unsigned int>& pos);
        std::vector<Node*> getLayerNodes(int layerIndex);
        int getMaxWeightedEdge(const std::string& nodeID, const std::string& edgeColor, const std::string& direction);
        uint getEdgeWeightBetweenNodes(const std::string& node1, const std::string& node2, const std::string& edgeColor="");
        uint longEdgeSupport(const std::vector<std::string>& path, const Node* candidateID);
        std::string idxLayer2NodeID(uint layerIdx, uint nodeIdx);
        
        void populateGraph(std::string bcf_file_path);
        void connectNodeGroupAdj(const std::vector<Node *>& nodes, std::string edgeColor);
        void connectNodeGroupAll(const std::vector<Node *>& nodes, std::string edgeColor, std::string readName);

        void printNodes();
        void printAdjList();
        void printConnectivityStats();
        void exportToDot();
        void exportToCytoscapeJSON();
        
        std::map<std::string, Node *> getNodes();
        void simulatePath(std::vector<std::string>& path, std::vector<uint>& choices, uint layer=0);
    private:
        static unsigned int nodeCount;
        unsigned int k;
        std::vector<std::vector<std::string>> layers;
        std::map<std::string, Node *> nodes;
        std::map<
            std::string, 
            std::pair<
                std::map< // stores out edges (connected to the SNVs in greater reference positions)
                    std::string, // name of the neighbor node
                    std::map<std::string, int> // color and weight
                >, 
                std::map< // stores in edges
                    std::string,
                    std::map<std::string, int>
                >
            >
        > adjList;

        std::vector<std::string> orderedNodes;
        std::vector<uint> orderedSNVs;
        std::map<std::string, std::vector<std::string>> readSupports; // map of read name to list of node IDs it supports
};

#endif