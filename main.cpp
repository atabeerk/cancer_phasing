#include <iostream>
#include <vector> 
#include <algorithm>
#include <string>

#include "Node.hpp"
#include "Graph.hpp"
#include "connectSNVs.hpp"
#include "util.hpp"


int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cout << "Usage: " << argv[0] << " <VCF file> <list of BAM files> <output file>" << std::endl;
        return 1;
    }

    // 1. VCF file (first argument)
    std::string vcfFile = argv[1];
    std::cout << "VCF File: " << vcfFile << std::endl;

    // 2. List of BAM files (second argument onwards, excluding last argument)
    std::vector<std::string> bamFiles;
    for (int i = 2; i < argc - 1; ++i) {
        bamFiles.push_back(argv[i]);
    }

    // Print BAM files
    if (bamFiles.size() > 6) {
        std::cout << "WARNING! only " << edgeColors.size() << " colors"
                  << "are used to draw graph edges but you provided "
                  << bamFiles.size() << " bam files." << std::endl;
    }
    std::cout << "BAM Files:" << std::endl;
    for (const auto& bamFile : bamFiles) {
        std::cout << "  - " << bamFile << std::endl;
    }

    // 3. Output file (last argument)
    std::string outputFile = argv[argc - 1];
    std::cout << "Output File: " << outputFile << std::endl;

    
    Graph g = Graph(1);
    g.populateGraph(vcfFile);
    for (uint i = 0; i < bamFiles.size(); i++){
        findConnectedSNVs(bamFiles[i], g, edgeColors[i % edgeColors.size()]);
    }

    g.exportToDot(outputFile);

    return 0;
}