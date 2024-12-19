#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <algorithm>

#include "util.hpp"

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
