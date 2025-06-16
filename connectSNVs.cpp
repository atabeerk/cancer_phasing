#include <iostream>
#include <vector> 
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <set>
#include <cstdio>

#include "htslib/sam.h"

#include "connectSNVs.hpp"
#include "Node.hpp"
#include "Graph.hpp"
#include "util.hpp"
#include "config.hpp"


bool readSupportsAllele(const ReadInfo& read, const Node& node) {
    bool support = false;

    for (const auto& snv : read.SNVs) {
        for (size_t i = 0; i < node.posArr().size(); i++) {
            if (snv.pos == std::stoi(node.posArr()[i])) { // this can be improved
                char read_base = std::toupper(static_cast<unsigned char>(snv.base[0]));
                char allele_base = std::toupper(static_cast<unsigned char>(node.baseArr()[i][0]));
                if (read_base == allele_base) {
                    support = true;
                } else {
                    return false;
                }
            }
        }
    }
    return support;
}


void connectGraphNodes(const std::map<std::string, ReadInfo>& read_data, Graph& g, std::string edgeColor) {
    size_t i = 0;
    for (auto it = read_data.begin(); it != read_data.end(); ++it) {
        // Get the key (read name) and value (ReadInfo) from the iterator
        const std::string& readName = it->first;
        const ReadInfo& read = it->second;

        int k = g.getK();
        std::vector<Node *> nodesInRange;
        std::vector<uint> SNVsInRange;
        for (const auto& snv : read.SNVs) {
            SNVsInRange.push_back(snv.pos);
        }

        std::vector<uint> k_SNVs;
        for (int i = 0; i < SNVsInRange.size() - k + 1 && SNVsInRange.size() > k; i++) {
            k_SNVs.assign(SNVsInRange.begin() + i, SNVsInRange.begin() + i + k);

            std::vector<Node *> SNVnodes;
            // Get the two nodes corresponding to this SNV position
            try {
                SNVnodes = g.getSNVNodes(read.chrom, k_SNVs);
            } catch (const std::exception& e) {
                Config::getInstance().log("Exception while getting node for read: " + readName);
                Config::getInstance().log("Chromosome: " + read.chrom);
                Config::getInstance().log("Exception message: " + std::string(e.what()));
                Config::getInstance().log("SNVs: ");
                for (const auto& snv : k_SNVs) {
                    Config::getInstance().log("\tPosition: " + std::to_string(snv));
                }
            }

            nodesInRange.insert(nodesInRange.begin(), SNVnodes.begin(), SNVnodes.end());
        }
        std::vector<Node *> nodesToConnect;
        for (const auto& node : nodesInRange) {
            if (readSupportsAllele(read, *node)) {
                nodesToConnect.push_back(node);
            }
        }
        g.connectNodeGroup(nodesToConnect, edgeColor);
        i++;
    }
    std::cout << "Processed " << i << " reads for connection" << std::endl;
}


// Function to run pileup at SNV positions from the Graph's nodes
void findConnectedSNVs(const std::string& bamFilePath, Graph& graph, std::string edgeColor) {
    const auto& node_map = graph.getNodes();  // Access nodes using getNodes()
    std::vector<Node*> node_list;

    for (auto it = node_map.begin(); it != node_map.end(); ++it) {
        Node* node = it->second;  // Value of type Node*
        node_list.push_back(node);
    }

    // Iterate over SNV positions from nodes in the graph
    std::vector<std::string> regionLines;
    std::set<std::string> lineSet; // for keeping track of duplicates in regionLines
    for (const auto& node : node_list) {
        for (const std::string& pos_str : node->posArr()) {
            uint pos = std::stoi(pos_str) + 1; // convert to 1-based
             std::string region = node->chrom() + "\t" + std::to_string(pos);
             if (lineSet.find(region) == lineSet.end()) {
                lineSet.insert(region);
                regionLines.push_back(region);
             }
        }
    }

    std::string regionListFile = Config::getInstance().getOutputDir() + "/regionList.txt";
    std::string mpileupPath = Config::getInstance().getOutputDir() + "/mpileup.out";
    writeToFile(regionListFile, regionLines);
    executeMpileupCmd(bamFilePath, regionListFile, mpileupPath);
    std::map<std::string, ReadInfo> read_data = readMpileupOutput(mpileupPath);
    // printParsedMpileupOutput(read_data);
    connectGraphNodes(read_data, graph, edgeColor);
}


void executeMpileupCmd(std::string bamFilePath, std::string regionListFile, std::string mpileupPath) {
    // Run the command
    std::string command = "samtools mpileup -Q 0 -q 0 -O --output-QNAME --output-extra POS,RLEN -l " + regionListFile +
                          " " + bamFilePath + " -o " + mpileupPath;
    int ret_code = system(command.c_str());

    // Check the return code
    if (ret_code == 0) {
        std::cout << "mpileup command executed successfully." << std::endl;
    } else {
        std::cerr << "mpileup command execution failed with code: " << ret_code << std::endl;
    }

}


std::map<std::string, ReadInfo> readMpileupOutput(std::string filename) {
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Failed to open mpileup output the file!" << std::endl;
    }

    std::map<std::string, ReadInfo> read_data;  // Data structure to store read information

    std::string line;
    size_t j = 0;
    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::string chromosome, reference_position_str, bases, base_positions_str, read_names, skip, start_pos_str, rlen_str;

        // Read columns and skip some
        std::getline(ss, chromosome, '\t');             // Column 1: Chromosome
        std::getline(ss, reference_position_str, '\t'); // Column 2: Reference position
        std::getline(ss, skip, '\t');                   // Column 3: Reference base
        std::getline(ss, skip, '\t');                   // Column 4: Number of elements (skip this one)
        std::getline(ss, bases, '\t');                  // Column 5: Bases of the reads
        std::getline(ss, skip, '\t');                   // Column 6: Base quality
        std::getline(ss, skip, '\t');                   // Column 7: Read position
        std::getline(ss, read_names, '\t');             // Column 8: Read names (comma-separated)
        std::getline(ss, start_pos_str, '\t');          // Column 9: start positions (comma-separated)
        std::getline(ss, rlen_str, '\t');               // Column 10: read lengths (comma-separated)


        // Convert reference position to integer and make 0-based
        int reference_position = std::stoi(reference_position_str) - 1;

        // Split read_names by commas
        std::istringstream read_names_stream(read_names);
        std::istringstream start_pos_stream(start_pos_str);
        std::istringstream rlen_stream(rlen_str);
        std::string read_name, start_pos, rlen;
        size_t i = 0;
        while (std::getline(read_names_stream, read_name, ',') && 
               std::getline(start_pos_stream, start_pos, ',') && 
               std::getline(rlen_stream, rlen, ',')) {
            SnvInfo snv_info;
            snv_info.pos = reference_position;
            snv_info.base = bases[i];
            if (read_data.find(read_name) != read_data.end()) {
                read_data[read_name].SNVs.push_back(snv_info);
            } else{
                ReadInfo read_info;
                try {
                    read_info.chrom = chromosome;
                    read_info.start_pos = std::stoi(start_pos);
                    read_info.rlen = std::stoi(rlen);
                    read_info.SNVs.push_back(snv_info);
                    read_data[read_name] = read_info;
                } catch (...) {
                    std::cerr << "Exception at read: " << i << std::endl;
                }
            }
            i++;
        }
        j++;
    }

    // Close the file after processing
    file.close();
    return read_data;
}


void printParsedMpileupOutput(const std::map<std::string, ReadInfo>& read_data) {
    size_t c=0;
    for (const auto& outer_entry : read_data) {
        std::cout << "Read name: " << outer_entry.first << std::endl;  // Outer map key (read name)
        
        const ReadInfo& read_info = outer_entry.second;  // ReadInfo struct

        // Print the members of ReadInfo
        std::cout <<  "Chromosome: " << read_info.chrom
                  << "  Start Position: " << read_info.start_pos << ", Read Length: " << read_info.rlen << std::endl;

        // Iterate over the SNVs for this read
        for (const auto& snv : read_info.SNVs) {
            c++;
            // Print each SnvInfo in a single line
            std::cout << "  SNV: " 
                      << "Position: " << snv.pos << ", "
                      << "Base: " << snv.base << std::endl;
        }
    }
    std::cout << "Total number of entries (positions X reads per position): " << c << std::endl;
}
