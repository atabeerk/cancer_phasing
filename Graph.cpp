#include <iostream>
#include <fstream>
#include <algorithm> 

#include "Graph.hpp"
#include "util.hpp"

unsigned int Graph::nodeCount;

Graph::Graph(unsigned int k) {
    /* 
    k is the number of SNVs each node will store
    */
    this->k = k;
    this->nodeCount = 0;
}


void Graph::addSNV(std::string chrom, std::vector<unsigned int> pos, std::vector<std::string> ref, std::vector<std::string> alt) {
    /*
    Adding an SNV to the graph constitutes of adding two independent nodes,
    one for each hapltype (ref/alt). Positioins of the two node objects will be
    the same but bases on these positions will be different.
    */

    addNode(chrom, pos, ref, "H1");
    addNode(chrom, pos, alt, "H2");

    orderedSNVs.push_back(pos[0]);
}


void Graph::addNode(std::string chrom, std::vector<unsigned int> pos, std::vector<std::string> bases, std::string suffix) {
    Node *n = new Node(k, chrom, pos, bases, suffix);
    nodes[n->ID()] = n;

    // Create an empty entry in the adjaceny list
    adjList[n->ID()] = std::map<std::string, int>();
    orderedNodes.push_back(n->ID());

    nodeCount++;
}


std::vector<Node*> Graph::getSNVNodes(std::string chrom, std::vector<unsigned int> pos){
    /*
    Returns both nodes that correspond to the position in the nodeID
    */
    Node* node1 = getNode(chrom, pos, "H1");
    Node* node2 = getNode(chrom, pos, "H2");
    std::vector<Node*> v = {node1, node2};
    return v;
}


Node* Graph::getNode(std::string chrom, std::vector<unsigned int> pos, std::string suffix){
    std::string nodeID = createNodeID(chrom, pos, suffix);
    return nodes[nodeID];
}


Node* Graph::getNode(std::string nodeID){
    return nodes[nodeID];
}


unsigned int Graph::getK(){
    return k;
}


unsigned int Graph::getNodeCount(){
    return nodeCount;
}


void Graph::addEdge(std::string node1, std::string node2, int weight){
    if (adjList[node1].find(node2) != adjList[node1].end()) {
        adjList[node1][node2] += weight;    
    }
    else {
        adjList[node1][node2] = weight; 
    }
}


void Graph::connectNodeGroup(const std::vector<Node *>& nodes){
    /*
    Only connects subsequentnodes in the list
    E.g., For [n1, n2, n3],  n1-n2 and n2-n3 are connencted
    */

    if (nodes.size() < 2) {
        return;
    }

    for(int i = 0; i < nodes.size() - 1; i++) {
            addEdge(nodes[i]->ID(), nodes[i+1]->ID());
    }
}


void Graph::printNodes(){
    std::cout << nodeCount << " nodes in total" <<std::endl;
    for (const auto& pair : nodes){
        pair.second->print();
    }
}


void Graph::printAdjList() {
    // Iterate through the outer map
    for (const auto& entry : adjList) {
        std::string node = entry.first;  // Outer map key (the node)
        const std::map<std::string, int>& neighbors = entry.second;  // Inner map (neighbors and weights)

        std::cout << "Node " << node << " has edges to:\n";

        // Iterate through the inner map (neighbors and weights)
        for (const auto& neighbor : neighbors) {
            std::cout << "  Neighbor: " << neighbor.first << ", Weight: " << neighbor.second << "\n";
        }
    }
}


void Graph::exportToDot(const std::string& filename) {
    std::ofstream file(filename);
    file << "graph G {\n";  // use "graph G" for undirected, "digraph G" for directed
    file << "    rankdir=LR;\n";

    int h1_index = 0, h2_index = 0;  // Track x-coordinates for H1 and H2 nodes
    const int h1_y = 1;  // Fixed y-coordinate for H1 line
    const int h2_y = -1; // Fixed y-coordinate for H2 line

    // Assign positions and colors to H1 and H2 nodes
    for (const auto& node : orderedNodes) {
        if (endsWith(node, "H1")) {
            file << "    \"" << node << "\" [pos=\"" << h1_index * 3 << "," << h1_y << "!\", color=blue, style=filled, fillcolor=lightblue];\n";
            h1_index++;  // Increment x-coordinate for the next H1 node
        } else if (endsWith(node, "H2")) {
            file << "    \"" << node << "\" [pos=\"" << h2_index * 3 << "," << h2_y << "!\", color=red, style=filled, fillcolor=pink];\n";
            h2_index++;  // Increment x-coordinate for the next H2 node
        }
    }
    
    // Iterate over the adjacency list
    for (std::vector<std::string>::iterator it = orderedNodes.begin(); it != orderedNodes.end(); ++it) {
        const std::string& nodeID = *it;
        std::string nodePos = getNode(nodeID)->posArr()[0];
        const std::map<std::string, int>& neighbors = adjList[nodeID];

        // Iterate through neighbors
        for (std::map<std::string, int>::const_iterator jt = neighbors.begin(); jt != neighbors.end(); ++jt) {
            const std::string& neighborID = jt->first;
            std::string neighborPos = getNode(neighborID)->posArr()[0];

            // Only create edges between subsequent SNV positions. Because there are
            // H1 and H2 nodes, we check if difference between the indices is less than 2 (instead of 1)
            if (distanceBetweenElements(orderedSNVs, nodePos, neighborPos) < 2) {
                int weight = jt->second;
                // Writing edges with weight as label
                file << "    \"" << nodeID << "\" -- \"" << neighborID << "\" [label=\"" << weight << "\"];\n";
            }
        }
    }
    file << "}\n";
    file.close();

    std::cout << "Graph exported to " << filename << std::endl;
}

