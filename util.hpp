#ifndef UTIL_H
#define UTIL_H

#include <string>
#include <vector>
#include <iostream>

extern std::vector<std::string> edgeColors; // used to color edges in .dot
struct Cluster {
    std::vector<size_t> indices; // indices of paths in this cluster
};

std::string joinStringVector(const std::vector<std::string>& elements, const std::string& separator="-");
std::vector<std::string> uintToStringVector(std::vector<uint> s);
std::string createNodeID(std::string chrmo, std::vector<uint> pos, std::string suffix="");
bool endsWith(const std::string& str, const std::string& suffix);
int findIndex(std::vector<std::string>vec , std::string item);
int distanceBetweenElements(const std::vector<uint>& list, uint elem1, uint elem2);
int distanceBetweenElements(const std::vector<uint>& list, std::string elem1, std::string elem2);
void writeToFile(std::string filename, std::vector<std::string> content);
std::vector<double> softmax_with_temperature(const std::vector<double>& logits, double T=1.0);
double hammingDistance(const std::vector<uint>& a, const std::vector<uint>& b);
double clusterDistance(const Cluster& c1, const Cluster& c2, const std::vector<std::vector<double>>& D);
std::vector<Cluster> hierarchicalClustering(const std::vector<std::vector<uint>>& paths, double threshold);
#endif