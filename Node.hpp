#ifndef NODE_H
#define NODE_H

#include <string>
#include <vector>


class Node {
    public:
        Node(char k, std::string chrom, std::vector<unsigned int> pos, std::vector<std::string> bases, std::string suffix="");
        std::string ID();
        const std::vector<std::string>& posArr();
        const std::vector<std::string>& baseArr();
        const std::string& chrom();
        ~Node();
        void print();
    private:
        char k;
        std::string nodeID;
        std::string chr;
        std::vector<std::string> pos_arr;
        std::vector<std::string> base_arr;

};

#endif