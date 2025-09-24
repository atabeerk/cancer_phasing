#include <iostream>
#include <vector> 
#include <algorithm>
#include <string>

#include <cxxopts.hpp>

#include "Node.hpp"
#include "Graph.hpp"
#include "connectSNVs.hpp"
#include "util.hpp"
#include "config.hpp"



int main(int argc, char* argv[]) {
    std::string vcfPath;
    std::vector<std::string> bamFiles;
    std::string outDir;
    int k;
    try {
        cxxopts::Options options("CancerPhaser", "de Bruijn inspired phasing of cancer clones using SNVs and BAM files");

        options.add_options()
            ("vcf", "Path to VCF file", cxxopts::value<std::string>())
            ("bam", "Paths to BAM files", cxxopts::value<std::vector<std::string>>())
            ("out-dir", "Output directory", cxxopts::value<std::string>())
            ("k", "K value", cxxopts::value<int>())
            ("h,help", "Print help");

        auto result = options.parse(argc, argv);

        if (result.count("help")) {
            std::cout << options.help() << std::endl;
            return 0;
        }

        vcfPath = result["vcf"].as<std::string>();
        bamFiles = result["bam"].as<std::vector<std::string>>();
        outDir = result["out-dir"].as<std::string>();
        k = result["k"].as<int>();

    } catch (const std::exception& e) {
        std::cerr << "Error parsing options: " << e.what() << std::endl;
        return 1;
    }

    Config::getInstance().setOutputDir(outDir);

    // Log input parameters
    Config::getInstance().log("=== Program Inputs ===\n");
    Config::getInstance().log("VCF: " + vcfPath + "\n");
    Config::getInstance().log("Output Directory: " + outDir + "\n");
    Config::getInstance().log("k: " + std::to_string(k) + "\n");
    Config::getInstance().log("BAM files:\n");
    for (const auto& bam : bamFiles) {
        Config::getInstance().log("  " + bam + "\n");
    }
    Config::getInstance().log("======================\n");

    // Start processing
    Graph g = Graph(k);
    g.populateGraph(vcfPath);
    for (uint i = 0; i < bamFiles.size(); i++){
        findConnectedSNVs(bamFiles[i], g, edgeColors[i % edgeColors.size()]);
    }

    g.printConnectivityStats();
    g.exportToCytoscapeJSON();
    g.exportToDot();

    return 0;
}