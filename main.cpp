#include <iostream>
#include <vector> 
#include <algorithm>
#include <string>

#include "Node.hpp"
#include "Graph.hpp"
#include "connectSNVs.hpp"
#include "util.hpp"
#include "config.hpp"


int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cout << "Usage: " << argv[0] << " <VCF file> <list of BAM files> <output folder>" << std::endl;
        return 1;
    }

    // 1. VCF file (first argument)
    std::string vcfFile = argv[1];
    std::cout << "VCF File:\n  - " << vcfFile << std::endl;

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

    // 3. Output folder (last argument)
    std::string dir = argv[argc - 1];
    std::cout << "Output Folder: " << dir << std::endl;

    Config::getInstance().setOutputDir(dir);

    Graph g = Graph(5);
    g.populateGraph(vcfFile);
    for (uint i = 0; i < bamFiles.size(); i++){
        findConnectedSNVs(bamFiles[i], g, edgeColors[i % edgeColors.size()]);
    }

    g.exportToDot();

    return 0;
}