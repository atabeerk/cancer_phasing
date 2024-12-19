#include <iostream>

#include "Node.hpp"
#include "util.hpp"


Node::Node(char k, std::string chrom, std::vector<unsigned int> pos, std::vector<std::string> bases, std::string suffix) {
    for(int i=0; i < k; i++){
        this->pos_arr.push_back(std::to_string(pos[i]));
        this->base_arr.push_back(bases[i]);
    }
    this->k = k;
    std::string nodeID = createNodeID(chrom, pos, suffix);
    this->nodeID = nodeID;
}


std::string Node::ID(){
    return nodeID;
}


const std::vector<std::string>& Node::posArr(){
    return pos_arr;
}


const std::vector<std::string>& Node::baseArr(){
    return base_arr;
}


const std::string& Node::chrom() {
    return chr;
}


void Node::print() {
    for (auto i : pos_arr){
        std::cout << i << " ";
    }
    std::cout << std::endl;

    for (auto i : base_arr){
        std::cout << i << " ";
    }
    std::cout << std::endl;

}

Node::~Node(){

}