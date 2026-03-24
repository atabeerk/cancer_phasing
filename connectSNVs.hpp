#ifndef SNVS_H
#define SNVS_H

#include "Graph.hpp"

struct SnvInfo {
    uint pos;
    std::string base;
};

struct ReadInfo {
    std::string chrom;
    uint start_pos;
    uint rlen;
    std::vector<SnvInfo> SNVs;
};

bool readSupportsAllele(const ReadInfo& read, const Node& node);
void connectGraphNodes(const std::map<std::string, ReadInfo>& read_data, Graph& g, std::string edgeColor);
void findConnectedSNVs(const std::string& bam_path, Graph& graph, std::string edgeColor="black");
void executeMpileupCmd(std::string bamFilePath, std::string regionListFile, std::string mpileupPath);
std::map<std::string, ReadInfo> readMpileupOutput(std::string filename);
void printParsedMpileupOutput(const std::map<std::string, ReadInfo>& read_data);
#endif