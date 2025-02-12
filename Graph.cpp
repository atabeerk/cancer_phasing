#include <iostream>
#include <fstream>
#include <algorithm> 
#include <cstdio>

#include "htslib/vcf.h"
#include "htslib/sam.h"

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
    std::pair<std::map<std::string,std::map<std::string, int>>, std::map<std::string,std::map<std::string, int>>> emptyPair;
    adjList[n->ID()] = emptyPair;
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
    return nodes.at(nodeID);
    
}


Node* Graph::getNode(std::string nodeID){
    return nodes[nodeID];
}


std::vector<uint> Graph::getOrderedSNVs() {
    return orderedSNVs;
}

// TODO: store this value for each edge for efficiency
int Graph::getMaxWeightedEdge(const std::string& nodeID, const std::string& edgeColor, const std::string& direction) {
    const auto& neighbors = (direction == "out")
        ? adjList.at(nodeID).first
        : adjList.at(nodeID).second;

    int maxValue = -1;
    bool found = false;

    for (const auto& n : neighbors) {  // Iterate over the outer map
        auto edges = n.second;
        auto it = edges.find(edgeColor);  // Check if edgeColor exists in the inner map
        if (it != edges.end()) {
            maxValue = std::max(maxValue, it->second);  // Track the max value
            found = true;
        }
    }

    if (!found) {
        throw std::runtime_error("Error: Edge color not found in any inner map.");
    }

    return maxValue;
}


/*
Given a vcf file, populates the graph with nodes corresponding to SNVs.
Each SNV position corresponds to two nodes (one for each haplotype).
*/
void Graph::populateGraph(std::string bcf_file_path){
    htsFile *test_bcf = NULL;
    bcf_hdr_t *test_header = NULL;
    bcf1_t *test_record = bcf_init();
    test_bcf = bcf_open(bcf_file_path.c_str(), "r");
    if(test_bcf == NULL) {
        throw std::runtime_error("Unable to open file.");
    }
    test_header = bcf_hdr_read(test_bcf);
    if(test_header == NULL) {
        throw std::runtime_error("Unable to read header.");
    }

    uint j = 0;
    std::vector<unsigned int> pos_arr;
    std::vector<std::string> ref_arr;
    std::vector<std::string> alt_arr;

    unsigned int k = getK();
    while(bcf_read(test_bcf, test_header, test_record) == 0) {
        j++;
        bcf_unpack(test_record, BCF_UN_ALL);

        if (j % 100000 == 0){
            std::cout << "Processed " << j << " records" << std::endl;
        }

        const char *filter_string = bcf_hdr_int2id(test_header, BCF_DT_ID, test_record->d.flt[0]);
        if (strcmp(filter_string, "PASS") != 0) {
            continue; // Skip records that do not have "PASS" in the FILTER column
        }

        // Extract REF and ALT alleles
        std::string ref = test_record->d.allele[0];
        std::string alt = test_record->d.allele[1];

        // Check if the variant is a SNP (both REF and ALT must have length 1)
        if (ref.length() != 1 || alt.length() != 1) {
            continue; // Skip non-SNPs
        }

        std::string region = bcf_hdr_id2name(test_header, test_record->rid);

        pos_arr.push_back(test_record->pos);
        ref_arr.push_back(ref);
        alt_arr.push_back(alt);

        // For the first k-1 SNVs, there are not enough # of SNVs to create a node
        if (j < k) {
            continue;
        }

        addSNV(region, pos_arr, ref_arr, alt_arr);

        // Slide the window by removing the first SNV in the arrays
        pos_arr.erase(pos_arr.begin());
        ref_arr.erase(ref_arr.begin());
        alt_arr.erase(alt_arr.begin());
    }

    bcf_hdr_destroy(test_header);
    bcf_destroy(test_record); 
    bcf_close(test_bcf);
    std::cout << ".vcf file parsed, number of nodes in the graph: ";
    std::cout << getNodeCount() << std::endl;
}


unsigned int Graph::getK(){
    return k;
}


unsigned int Graph::getNodeCount(){
    return nodeCount;
}


std::map<std::string, Node *> Graph::getNodes(){
    return nodes;
}


void Graph::addEdge(std::string node1, std::string node2, int weight, std::string edgeColor, std::string direction){
    std::map<std::string, std::map<std::string, int>>& pairElem = (direction == "out")
    ? adjList.at(node1).first
    : adjList.at(node1).second;
    if (pairElem.count(node2) > 0) {

        if (pairElem[node2].count(edgeColor) > 0) {
            pairElem[node2][edgeColor] += weight;    
        }
        else {
            pairElem[node2][edgeColor] = weight;
        }
    }
    else {
        pairElem[node2] = std::map<std::string, int>();
        pairElem[node2][edgeColor] = weight;
    }
}


void Graph::connectNodeGroup(const std::vector<Node *>& nodes, std::string edgeColor){
    /*
    Only connects subsequentnodes in the list
    E.g., For [n1, n2, n3],  n1-n2 and n2-n3 are connencted
    */

    if (nodes.size() < 2) {
        return;
    }

    for(int i = 0; i < nodes.size() - 1; i++) {
            addEdge(nodes[i]->ID(), nodes[i+1]->ID(), 1, edgeColor, "out");
            addEdge(nodes[i+1]->ID(), nodes[i]->ID(), 1, edgeColor, "in");
    }
}


void Graph::printNodes(){
    std::cout << nodeCount << " nodes in total" <<std::endl;
    for (const auto& pair : nodes){
        pair.second->print();
    }
}


void Graph::printAdjList() {
    // std::cout << "Graph Representation (Adjacency List):\n";

    // for (const auto& node : adjList) {
    //     std::cout << "Node: " << node.first << "\n";
        
    //     for (const auto& neighbor : node.second) {
    //         std::cout << "  -> " << neighbor.first << " (";
            
    //         bool first = true;
    //         for (const auto& edge : neighbor.second) {
    //             if (!first) std::cout << ", ";
    //             std::cout << "Color: " << edge.first << ", Weight: " << edge.second;
    //             first = false;
    //         }
    //         std::cout << ")\n";
    //     }
    //     std::cout << "---------------------------------\n";
    // }
}


void Graph::exportToDot(const std::string& filename) {
    std::string tmp_file = "tmp_file.dot";
    std::ofstream file(tmp_file);
    file << "graph G {\n";  // use "graph G" for undirected, "digraph G" for directed
    // file << "   graph [dpi = 300];\n";
    file << "   rankdir=LR;\n";

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
        const std::map<std::string, std::map<std::string, int>> outNeighbors = adjList.at(nodeID).first;

        // Iterate through neighbors
        for (std::map<std::string, std::map<std::string, int>>::const_iterator jt = outNeighbors.begin(); jt != outNeighbors.end(); ++jt) {
            
            const std::string& neighborID = jt->first;
            std::string neighborPos = getNode(neighborID)->posArr()[0];

            const std::map<std::string, int> edges = jt->second;

            // Only create edges between subsequent SNV positions. Because there are
            // H1 and H2 nodes, we check if difference between the indices is less than 2 (instead of 1)
            if (distanceBetweenElements(orderedSNVs, nodePos, neighborPos) < 2) {
                // Iterate over out edges from the current node, jt-> first stores the out edges.
                for (std::map<std::string, int>::const_iterator  kt = edges.begin(); kt != edges.end(); ++kt) {
                    std::string edgeColor = kt->first;
                    int weight = kt->second;

                    int maxOutWeight = getMaxWeightedEdge(nodeID, edgeColor, "out");
                    int maxInWeight = getMaxWeightedEdge(neighborID, edgeColor, "in");
                    
                    if (weight * 5 >= maxOutWeight) {
                        // Writing edges with weight as label
                        file << "    \"" << nodeID << "\" -- \"" << neighborID << "\" [label=\"" << weight << "\" color=\"" << edgeColor << "\"];\n";
                    }
                }
            }
        }
    }
    file << "}\n";
    file.close();

    
    std::string png_convert_cmd = "dot -Tsvg " + tmp_file + " -o " + filename;
    if (system(png_convert_cmd.c_str()) != 0) {
        std::cerr << "Error converting graph from .dot" << std::endl;
    } else {
        std::cout << "Graph exported to " << filename << std::endl;
    }
    // std::remove(tmp_file.c_str());

    
}

