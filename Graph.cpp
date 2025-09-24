#include <iostream>
#include <fstream>
#include <algorithm> 
#include <cstdio>

#include "htslib/vcf.h"
#include "htslib/sam.h"

#include "json.hpp"
using json_t = nlohmann::json;

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


void logHaplotypeMap(const std::vector<std::string>& bases,
                     const std::vector<std::string>& pos,
                     const std::string& hapID) {
    Config::getInstance().log(hapID);
    for (size_t i = 0; i < bases.size(); ++i) {
        Config::getInstance().log("\t" + pos[i] + ":" + bases[i]);
    }
    Config::getInstance().log("\n");
}


void Graph::addSNV(const std::string& chrom,
                   const std::vector<unsigned int>& pos,
                   const std::vector<std::string>& ref,
                   const std::vector<std::string>& alt) {
    /*
    Adds all combinations of ref/alt alleles as separate nodes to the graph.
    Each combination represents a possible haplotype.
    */

    if (pos.empty() || ref.size() != pos.size() || alt.size() != pos.size()) {
        throw std::invalid_argument("SNV input vectors must be the same size and non-empty.");
    }

    size_t n = pos.size();
    size_t num_combinations = 1 << n; // 2^n combinations

    for (size_t i = 0; i < num_combinations; ++i) {
        std::vector<std::string> hap;
        for (size_t j = 0; j < n; ++j) {
            // If bit j is set, use alt[j]; else, use ref[j]
            if (i & (1 << j)) {
                hap.push_back(alt[j]);
            } else {
                hap.push_back(ref[j]);
            }
        }
        // Label each combination with a haplotype ID like H0, H1, ...
        std::string hapID = "H" + std::to_string(i);
        addNode(chrom, pos, hap, hapID);
    }

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
    logHaplotypeMap(n->baseArr(), n->posArr(), n->ID());
}


std::vector<Node*> Graph::getSNVNodes(const std::string& chrom, const std::vector<unsigned int>& pos) {
    std::vector<Node*> nodes;
    size_t num_haplotypes = 1 << this->k;  // 2^k

    for (size_t i = 0; i < num_haplotypes; ++i) {
        std::string hapID = "H" + std::to_string(i);
        Node* node = getNode(chrom, pos, hapID);
        if (node != nullptr) {
            nodes.push_back(node);
        }
    }

    return nodes;
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
Used in the Graph::populateGraph() function to determine if a record is
heterozygous and filter it out if not. This function checks the genotype (GT)
field of the first two samples in a BCF record. 
*/
bool any_sample_is_heterozygous(bcf_hdr_t* hdr, bcf1_t* rec) {
    int* gt_arr = nullptr;
    int ngt_arr = 0;

    // Get GT field
    int ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
    if (ngt <= 1) {
        free(gt_arr);
        return false;
    }

    int nsamples = bcf_hdr_nsamples(hdr);

    for (int i = 0; i < nsamples; ++i) {
        int allele1 = bcf_gt_allele(gt_arr[i * 2]);
        int allele2 = bcf_gt_allele(gt_arr[i * 2 + 1]);

        // Heterozygous: both are numbers and not equal
        if (allele1 >= 0 && allele2 >= 0 && allele1 != allele2) {
            free(gt_arr);
            return true;
        }
    }

    free(gt_arr);
    return false;
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
        throw std::runtime_error("Unable to open BCF file.");
    }
    test_header = bcf_hdr_read(test_bcf);
    if(test_header == NULL) {
        throw std::runtime_error("Unable to read BCF header.");
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

        if ( // Skip non-SNPs, low quality, or non-heterozygous
            ref.length() != 1 ||
            alt.length() != 1 ||
            test_record->qual < 20 
            // !any_sample_is_heterozygous(test_header, test_record)
        ) {
            continue;
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
    Config::getInstance().log("VCF file parsed, number of nodes in the graph: " + getNodeCount());
    Config::getInstance().log("Nodes:");
    for (const auto& pair : nodes) {
        Config::getInstance().log("\t" + pair.first);
    }
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

void Graph::connectNodeGroupAll(const std::vector<Node *>& nodes, std::string edgeColor){
    /*
    Connects all nodes in the list
    */
    if (nodes.empty()) {
        return;
    }

    for (size_t i = 0; i < nodes.size(); i++) {
        for (size_t j = i + 1; j < nodes.size(); j++) {
            addEdge(nodes[i]->ID(), nodes[j]->ID(), 1, edgeColor, "out");
            addEdge(nodes[j]->ID(), nodes[i]->ID(), 1, edgeColor, "in");
        }
    }
}


void Graph::connectNodeGroupAdj(const std::vector<Node *>& nodes, std::string edgeColor){
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


#include <unordered_set>

void Graph::printConnectivityStats() {
    size_t n = 1 << k;  // 2^k nodes per group
    size_t totalGroups = orderedNodes.size() / n;

    std::unordered_set<std::string> connectedNodeIDs;
    unsigned int connectedGroups = 0;
    uint totalConnectedLength = 0;

    // total genomic region
    uint totalRegionLength = 0;
    if (!orderedNodes.empty()) {
        Node* firstNode = nodes.at(orderedNodes.front());
        Node* lastNode  = nodes.at(orderedNodes.back());
        totalRegionLength = std::stoul(lastNode->posArr().back()) -
                            std::stoul(firstNode->posArr().front()) + 1;
    }

    size_t ng = 0;
    while (ng + 1 < totalGroups) {
        size_t startCurr = n * ng;
        size_t endCurr = startCurr + n;
        size_t startNext = n * (ng + 1);
        size_t endNext = startNext + n;

        uint chainStart = std::stoul(nodes.at(orderedNodes[startCurr])->posArr().front());
        uint chainEnd = std::stoul(nodes.at(orderedNodes[startCurr])->posArr().back());
        bool chainActive = false;

        size_t currentNg = ng;

        while (currentNg + 1 < totalGroups) {
            bool groupsConnected = false;
            bool currentGroupConnected = false;

            startCurr = n * currentNg;
            endCurr = startCurr + n;
            startNext = n * (currentNg + 1);
            endNext = startNext + n;

            // iterate all nodes in current group
            for (size_t i = startCurr; i < endCurr; i++) {
                const std::string& currID = orderedNodes[i];
                Node* currNode = nodes.at(currID);

                auto it = adjList.find(currID);
                if (it == adjList.end()) continue;

                const auto& inEdges = it->second.second; // use inedges now

                // check each neighbor
                for (const auto& neighborKV : inEdges) {
                    const std::string& neighborID = neighborKV.first;
                    const auto& colorMap = neighborKV.second;

                    // print all edges
                    for (const auto& colorKV : colorMap) {
                        const std::string& edgeColor = colorKV.first;
                        int weight = colorKV.second;
                        std::cout << "Checking in-edge: " << currID
                                  << " <- " << neighborID
                                  << " (color: " << edgeColor
                                  << ", weight: " << weight << ")" << std::endl;
                    }

                    // check if neighbor is in the next group
                    for (size_t j = startNext; j < endNext; j++) {
                        if (orderedNodes[j] == neighborID) {
                            groupsConnected = true;
                            currentGroupConnected = true;

                            connectedNodeIDs.insert(currID);

                            Node* lastNodeNextGroup = nodes.at(orderedNodes[endNext - 1]);
                            chainEnd = std::stoul(lastNodeNextGroup->posArr().back());

                            std::cout << "Group " << currentNg
                                      << " is connected to group " << (currentNg + 1)
                                      << std::endl;
                        }
                    }
                }
            }

            if (groupsConnected) {
                chainActive = true;
                if (currentGroupConnected) connectedGroups++;
                currentNg++;  // move to next group in chain
            } else {
                break; // no more consecutive connections
            }
        }

        if (chainActive) {
            std::cout << "Connected chain from " << chainStart
                      << " to " << chainEnd << std::endl;
            totalConnectedLength += (chainEnd - chainStart + 1);
        }

        ng = currentNg + 1; // move to next unprocessed group
    }

    std::cout << "Connected nodes: " << connectedNodeIDs.size() << "/" << nodes.size()
              << " (" << (100.0 * connectedNodeIDs.size() / nodes.size()) << "%)" << std::endl;

    std::cout << "Connected groups: " << connectedGroups << "/" << totalGroups
              << " (" << (100.0 * connectedGroups / totalGroups) << "%)" << std::endl;

    std::cout << "Connected region length: " << totalConnectedLength
              << " / Total region length: " << totalRegionLength
              << " (" << (100.0 * totalConnectedLength / totalRegionLength) << "%)" << std::endl;
}


void Graph::exportToDot() {
    std::string output_folder = Config::getInstance().getOutputDir();
    std::string tmp_file = output_folder + "/graph_tmp.dot";
    std::ofstream file(tmp_file);
    file << "graph G {\n";  // use "graph G" for undirected, "digraph G" for directed
    file << "   layout=neato;\n";
    // file << "   overlap=false;\n";  // optional, to avoid overlapping nodes
    file << "   splines=spline;\n";   // optional, smooth edges
    // file << "   graph [dpi = 300];\n";
    // file << "   rankdir=RL;\n";

    // Step 1: Group nodes by SNV position in order of appearance in orderedNodes
    std::map<std::string, std::vector<std::string>> pos_to_nodes;
    std::vector<std::string> sortedSNVs;  // keep track of order

    for (const auto& node : orderedNodes) {
        std::string pos = getNode(node)->posArr()[0];

        // If this is the first time seeing this SNV position, record order
        if (pos_to_nodes.find(pos) == pos_to_nodes.end()) {
            sortedSNVs.push_back(pos);
        }

        pos_to_nodes[pos].push_back(node);
    }

    // Step 2: Layout nodes horizontally by SNV, vertically by haplotype
    const int x_spacing = 300;
    const int y_spacing = 100;

    int x_index = 0;
    for (const auto& snv : sortedSNVs) {
        const auto& nodeList = pos_to_nodes[snv];

        for (size_t i = 0; i < nodeList.size(); ++i) {
            const std::string& node = nodeList[i];
            int x = x_index * x_spacing;
            int y = i * y_spacing;

            file << "    \"" << node << "\" [pos=\"" << x << "," << y << "!\", style=filled, fillcolor=lightblue];\n";
        }
        x_index++;
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

            // Only create edges between subsequent SNV positions.
            int visRange = Config::getInstance().getMaxHopMultiplier() * k;
            uint16_t nodeDist = distanceBetweenElements(orderedSNVs, nodePos, neighborPos);
            // if (nodeDist <= visRange) {
                // Iterate over out edges from the current node, jt->first stores the out edges.
                for (std::map<std::string, int>::const_iterator  kt = edges.begin(); kt != edges.end(); ++kt) {
                    std::string edgeColor = kt->first;
                    int weight = kt->second;

                    int maxOutWeight = getMaxWeightedEdge(nodeID, edgeColor, "out");
                    int maxInWeight = getMaxWeightedEdge(neighborID, edgeColor, "in");
                    
                    // if (weight * 5 >= maxOutWeight && weight > 1) {
                    //     // Writing edges with weight as label
                    //     file << "    \"" << nodeID << "\" -- \"" << neighborID << "\" [label=\"" << weight << "\" color=\"" << edgeColor << "\"];\n";
                    // }
                file << "    \"" << nodeID << "\" -- \"" << neighborID << "\" [label=\"" << weight << "\" color=\"" << edgeColor << "\"];\n";
                }
            // }
        }
    }
    file << "}\n";
    file.close();

    std::cout << "Converting graph from DOT to SVG format..." << std::endl;
    std::string filename = Config::getInstance().getOutputDir() + "/graph.svg";
    std::string png_convert_cmd = "neato -n2 -Tsvg " + tmp_file + " -o " + filename;
    if (system(png_convert_cmd.c_str()) != 0) {
        std::cerr << "Error converting graph from .dot" << std::endl;
    } else {
        std::cout << "Graph exported to " << filename << std::endl;
    }
}

void Graph::exportToCytoscapeJSON() {
    std::string output_folder = Config::getInstance().getOutputDir();
    std::string json_file = output_folder + "/graph.json";
    std::ofstream file(json_file);

    if (!file.is_open()) {
        std::cerr << "Failed to open file for writing: " << json_file << std::endl;
        return;
    }

    json_t graph_json;
    json_t nodes_json = json_t::array();
    json_t edges_json = json_t::array();

    // Step 1: Group nodes by SNV position
    std::map<uint, std::vector<std::string>> pos_to_nodes;
    for (const auto& node : orderedNodes) {
        uint pos = static_cast<uint>(std::stoul(getNode(node)->posArr()[0]));
        pos_to_nodes[pos].push_back(node);
    }

    const int k = getK();
    const int maxHop = Config::getInstance().getMaxHopMultiplier();
    const int x_spacing = 500;
    const int y_spacing = 200;

    std::map<std::string, std::pair<int, int>> node_positions;
    int x_index = 0;

    for (const auto& snv : orderedSNVs) {
        if (pos_to_nodes.find(snv) == pos_to_nodes.end()) {
            ++x_index;
            continue;
        }

        const auto& nodeList = pos_to_nodes[snv];
        for (size_t i = 0; i < nodeList.size(); ++i) {
            const std::string& node = nodeList[i];
            int x = x_index * x_spacing;
            int y = i * y_spacing;
            node_positions[node] = {x, y};

            nodes_json.push_back({
                {"data", {{"id", node}}},
                {"position", {{"x", x}, {"y", y}}}
            });
        }
        ++x_index;
    }

    for (const auto& nodeID : orderedNodes) {
        const std::string& nodePos = getNode(nodeID)->posArr()[0];
        const auto& outNeighbors = adjList.at(nodeID).first;

        for (const auto& neighborEntry : outNeighbors) {
            const std::string& neighborID = neighborEntry.first;
            const std::string& neighborPos = getNode(neighborID)->posArr()[0];
            const auto& edges = neighborEntry.second;

            uint16_t nodeDist = distanceBetweenElements(orderedSNVs, nodePos, neighborPos);
            int visRange = maxHop * k;
            // if (nodeDist <= visRange) {
                for (const auto& edgeEntry : edges) {
                    std::string edgeColor = edgeEntry.first;
                    int weight = edgeEntry.second;

                    int maxOutWeight = getMaxWeightedEdge(nodeID, edgeColor, "out");
                    // if (weight * 5 >= maxOutWeight && weight > 1) {
                        std::string edgeID = nodeID + "_" + neighborID;

                        edges_json.push_back({
                            {"data", {
                                {"id", edgeID},
                                {"source", nodeID},
                                {"target", neighborID},
                                {"weight", weight},
                                {"color", edgeColor}
                            }}
                        });
                    // }
                }
            // }
        }
    }

    graph_json["elements"]["nodes"] = nodes_json;
    graph_json["elements"]["edges"] = edges_json;

    file << graph_json.dump(4);  // Pretty print
    file.close();

    std::cout << "Graph exported to Cytoscape JSON at " << json_file << std::endl;
}
