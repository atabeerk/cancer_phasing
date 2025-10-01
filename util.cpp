#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <limits>

#include "util.hpp"
#include "config.hpp"

std::vector<std::string> edgeColors = {"red", "blue", "green", "orange", "yellow", "purple"}; 

std::string joinStringVector(const std::vector<std::string>& elements, const std::string& separator) {
    std::string result;

    for (size_t i = 0; i < elements.size(); ++i) {
        result += elements[i];
        if (i < elements.size() - 1) {
            result += separator;
        }
    }

    return result;
}


std::vector<std::string> uintToStringVector(std::vector<uint> s){
    std::vector<std::string> result;
    for (uint e : s){
        result.push_back(std::to_string(e));
    }

    return result;
}

std::string createNodeID(std::string chrom, std::vector<uint> pos, std::string suffix) {
    std::vector<std::string> pos_string = uintToStringVector(pos);
    pos_string.insert(pos_string.begin(), chrom);
    pos_string.insert(pos_string.end(), suffix);
    std::string nodeID = joinStringVector(pos_string, "-");

    return nodeID;
}

bool endsWith(const std::string& str, const std::string& suffix) {
    if (str.size() < suffix.size()) {
        return false;
    }
    return str.substr(str.size() - suffix.size()) == suffix;
}

int findIndex(std::vector<std::string>vec , std::string elem) {
    auto it = std::find(vec.begin(), vec.end(), elem);

    // Check if the string was found
    if (it != vec.end()) {
        // Calculate the index by subtracting the iterator from the beginning of the vector
        return std::distance(vec.begin(), it);
    }
    return -1;
}

int distanceBetweenElements(const std::vector<uint>& list, uint elem1, uint elem2) {
    //shoudl be based on k
    // Find the iterators pointing to the elements
    auto it1 = std::find(list.begin(), list.end(), elem1);
    auto it2 = std::find(list.begin(), list.end(), elem2);

    // Check if both elements are found
    if (it1 == list.end() || it2 == list.end()) {
        return -1; // Return -1 if either element is not found
    }

    // Calculate and return the distance
    return std::abs(std::distance(it1, it2));
}

int distanceBetweenElements(const std::vector<uint>& list, std::string elem1, std::string elem2) {
    uint e1 = std::stoul(elem1);  // Converts string to unsigned int
    uint e2 = std::stoul(elem2); 
    return distanceBetweenElements(list, e1, e2);
}


void writeToFile(std::string filename, std::vector<std::string> content) {
    // Open the file, create if it doesn't exist
    std::ofstream outFile(filename);

    if (outFile.is_open()) {
        for (const auto& line : content) { // Assuming `lines` is a collection of strings
            outFile << line << '\n';
        }
        outFile.close();
    }
}


std::vector<double> softmax_with_temperature(const std::vector<double>& logits, double T) {
    std::vector<double> expVals;
    expVals.reserve(logits.size());

    std::vector<double> logits_double(logits.begin(), logits.end());

    double maxLogit = *max_element(logits_double.begin(), logits_double.end()); // for numerical stability
    double sum = 0.0;

    for (double z : logits_double) {
        double e = std::exp((z - maxLogit) / T);
        expVals.push_back(e);
        sum += e;
    }

    for (double& val : expVals) {
        val /= sum;
    }

    return expVals;
}


// Compute normalized Hamming distance between two integer vectors
double hammingDistance(const std::vector<uint>& a, const std::vector<uint>& b) {
    if (a.size() != b.size()) throw std::invalid_argument("Paths must have the same length");
    size_t diff = 0;
    for (size_t i = 0; i < a.size(); i++) if (a[i] != b[i]) diff++;
    return static_cast<double>(diff) / a.size();
}


// Compute the average distance between two clusters
double clusterDistance(const Cluster& c1, const Cluster& c2, const std::vector<std::vector<double>>& D) {
    double sum = 0.0;
    for (size_t i : c1.indices)
        for (size_t j : c2.indices)
            sum += D[i][j];
    return sum / (c1.indices.size() * c2.indices.size());
}


// Hierarchical clustering with distance threshold
std::vector<Cluster> hierarchicalClustering(const std::vector<std::vector<uint>>& paths, double threshold) {
    size_t n = paths.size();
    
    // Initial clusters: one path per cluster
    std::vector<Cluster> clusters(n);
    for (size_t i = 0; i < n; i++) clusters[i].indices.push_back(i);

    // Compute pairwise Hamming distances
    std::vector<std::vector<double>> D(n, std::vector<double>(n, 0.0));
    for (size_t i = 0; i < n; i++)
        for (size_t j = i + 1; j < n; j++) {
            D[i][j] = D[j][i] = hammingDistance(paths[i], paths[j]);
        }

    bool merged = true;
    while (merged && clusters.size() > 1) {
        merged = false;
        size_t mergeA = 0, mergeB = 0;
        double minDist = std::numeric_limits<double>::max();

        // Find the closest pair of clusters
        for (size_t i = 0; i < clusters.size(); i++) {
            for (size_t j = i + 1; j < clusters.size(); j++) {
                double dist = clusterDistance(clusters[i], clusters[j], D);
                if (dist < minDist) {
                    minDist = dist;
                    mergeA = i;
                    mergeB = j;
                }
            }
        }

        // Merge if below threshold
        if (minDist <= threshold) {
            clusters[mergeA].indices.insert(
                clusters[mergeA].indices.end(),
                clusters[mergeB].indices.begin(),
                clusters[mergeB].indices.end()
            );
            clusters.erase(clusters.begin() + mergeB);
            merged = true;
        }
    }

    return clusters;
}

