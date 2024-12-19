#include <iostream>
#include <vector> 
#include <algorithm>

#include <string>

#include "htslib/vcf.h"
#include "htslib/sam.h"

#include "Node.hpp"
#include "Graph.hpp"
#include "util.hpp"


void parse_vcf_file(char *bcf_file_path, Graph& g, std::vector<uint>& SNV_vector){
    htsFile *test_bcf = NULL;
    bcf_hdr_t *test_header = NULL;
    bcf1_t *test_record = bcf_init();
    test_bcf = bcf_open(bcf_file_path, "r");
    if(test_bcf == NULL) {
        throw std::runtime_error("Unable to open file.");
    }
    test_header = bcf_hdr_read(test_bcf);
    if(test_header == NULL) {
        throw std::runtime_error("Unable to read header.");
    }

    uint j = 0;
    std::vector<unsigned int> pos_arr;
    std::vector<std::string> ref_arr;
    std::vector<std::string> alt_arr;

    unsigned int k = g.getK();
    while(bcf_read(test_bcf, test_header, test_record) == 0) {
        j++;
        bcf_unpack(test_record, BCF_UN_ALL);

        if (j % 100000 == 0){
            std::cout << "Processed " << j << " records" << std::endl;
        }

        const char *filter_string = bcf_hdr_int2id(test_header, BCF_DT_ID, test_record->d.flt[0]);
        if (strcmp(filter_string, "PASS") != 0) {
            continue; // Skip records that do not have "PASS" in the FILTER column
        }

        // Extract REF and ALT alleles
        std::string ref = test_record->d.allele[0];
        std::string alt = test_record->d.allele[1];

        // Check if the variant is a SNP (both REF and ALT must have length 1)
        if (ref.length() != 1 || alt.length() != 1) {
            continue; // Skip non-SNPs
        }

        std::string region = bcf_hdr_id2name(test_header, test_record->rid);
        SNV_vector.push_back(test_record->pos);

        pos_arr.push_back(test_record->pos);
        ref_arr.push_back(ref);
        alt_arr.push_back(alt);

        // For the first k-1 SNVs, there are not enough # of SNVs to create a node
        if (j < k) {
            continue;
        }

        g.addSNV(region, pos_arr, ref_arr, alt_arr);

        // Slide the window by removing the first SNV in the arrays
        pos_arr.erase(pos_arr.begin());
        ref_arr.erase(ref_arr.begin());
        alt_arr.erase(alt_arr.begin());
    }

    bcf_hdr_destroy(test_header);
    bcf_destroy(test_record); 
    bcf_close(test_bcf);
    std::cout << ".vcf file parsed, number of nodes in the graph: ";
    std::cout << g.getNodeCount() << std::endl;
}


std::vector<std::string> getAlignments(const std::string& bamFile, const std::string& region) {
    std::vector<std::string> alignedReads;

    // Open the BAM file
    samFile *in = sam_open(bamFile.c_str(), "r");
    if (in == nullptr) {
        std::cerr << "Failed to open BAM file: " << bamFile << std::endl;
        return alignedReads;  // Return empty vector on error
    }

    // Load the BAM header
    bam_hdr_t *header = sam_hdr_read(in);
    if (header == nullptr) {
        std::cerr << "Failed to read BAM header." << std::endl;
        sam_close(in);
        return alignedReads;
    }

    // Load the BAM index
    hts_idx_t *idx = sam_index_load(in, bamFile.c_str());
    if (idx == nullptr) {
        std::cerr << "Failed to load BAM index." << std::endl;
        bam_hdr_destroy(header);
        sam_close(in);
        return alignedReads;
    }

    // Set up the iterator for the specified region
    hts_itr_t *iter = sam_itr_querys(idx, header, region.c_str());
    if (iter == nullptr) {
        std::cerr << "Failed to set up iterator for region: " << region << std::endl;
        hts_idx_destroy(idx);
        bam_hdr_destroy(header);
        sam_close(in);
        return alignedReads;
    }

    // Initialize a BAM record
    bam1_t *aln = bam_init1();
    
    // Iterate over the alignments in the specified region
    while (sam_itr_next(in, iter, aln) >= 0) {
        // Get the read name and add it to the vector
        alignedReads.push_back(bam_get_qname(aln));
    }

    // Clean up
    bam_destroy1(aln);
    hts_itr_destroy(iter);
    hts_idx_destroy(idx);
    bam_hdr_destroy(header);
    sam_close(in);

    return alignedReads;
}


std::vector<uint> extractRange(const std::vector<uint>& sortedVec, uint lower, uint upper) {
    // Find the first element that is >= lower
    auto startIt = std::lower_bound(sortedVec.begin(), sortedVec.end(), lower);
    
    // Find the first element that is > upper
    auto endIt = std::upper_bound(sortedVec.begin(), sortedVec.end(), upper);

    // Extract the range
    return std::vector<uint>(startIt, endIt); // endIt is exclusive
}


// Function to convert BAM record sequence to a string
std::string getReadSequence(const bam1_t* record) {
    const uint8_t* seq = bam_get_seq(record);
    int len = record->core.l_qseq;
    std::string seqStr;
    seqStr.reserve(len);

    for (int i = 0; i < len; ++i) {
        // BAM nucleotide encoding (1 = A, 2 = C, 4 = G, 8 = T, others for ambiguous bases)
        seqStr += seq_nt16_str[bam_seqi(seq, i)];
    }

    return seqStr;
}


bool doesReadSupportAllele(bam1_t* read, int snv_position, const std::string& allele) {
    // Convert 1-based VCF position to 0-based BAM position
    std::string readName = std::string(bam_get_qname(read));

    // Get the read's alignment information
    const bam1_core_t* core = &read->core;

    // Calculate 0-based position of the SNV relative to the reference
    int ref_start = core->pos; // 0-based position where the read starts on the reference
    int ref_end = ref_start + core->l_qseq; // Approximation of where it ends

    // If SNV is outside the range covered by the read, it cannot support the allele
    if (snv_position < ref_start || snv_position >= ref_end) {
        return false;
    }

    // Decode the CIGAR string to map reference positions to read positions
    uint32_t* cigar = bam_get_cigar(read);
    int ref_pos = ref_start;  // Reference position pointer
    int read_pos = 0;         // Read position pointer

    for (int i = 0; i < core->n_cigar; ++i) {
        uint32_t op_len = bam_cigar_oplen(cigar[i]);
        char op = bam_cigar_opchr(cigar[i]);

        switch (op) {
            case 'M': // Match or mismatch
            case '=': // Match
            case 'X': // Mismatch
                if (snv_position >= ref_pos && snv_position < ref_pos + op_len) {
                    // Determine the read position corresponding to the SNV
                    int offset = snv_position - ref_pos;
                    int read_index = read_pos + offset;

                    // Get the base from the read's sequence
                    uint8_t* seq = bam_get_seq(read);
                    char base = seq_nt16_str[bam_seqi(seq, read_index)];

                    // Compare with the provided allele
                    return std::toupper(base) == std::toupper(allele[0]);
                }
                ref_pos += op_len;
                read_pos += op_len;
                break;

            case 'I': // Insertion (does not advance the reference position)
                read_pos += op_len;
                break;

            case 'D': // Deletion (advances the reference position)
            case 'N': // Skipped region (advances the reference position)
                ref_pos += op_len;
                break;

            case 'S': // Soft clipping (does not advance the reference position)
            case 'H': // Hard clipping (does not advance the reference position)
                break;

            default:
                throw std::runtime_error("Unexpected CIGAR operation");
        }
    }

    // If no match was found, the read does not support the allele
    return false;
}


void connectGraphNodes(const char* bamFilePath, std::vector<uint>& SNVs, Graph& g) {
    // Open BAM file
    samFile* inFile = sam_open(bamFilePath, "r");
    if (inFile == nullptr) {
        std::cerr << "Failed to open BAM file: " << bamFilePath << std::endl;
        return;
    }

    // Load BAM header
    bam_hdr_t* header = sam_hdr_read(inFile);
    if (header == nullptr) {
        std::cerr << "Failed to read BAM header" << std::endl;
        sam_close(inFile);
        return;
    }

    // Initialize a BAM record
    bam1_t* record = bam_init1();

    // Iterate over each read in the BAM file
    uint matchedNodes = 0;
    uint readNo = 0;
    uint readXnodes = 0;
    while (sam_read1(inFile, header, record) >= 0) {
        readNo++;
        if (readNo % 1000 == 0) {
            std::cout << readNo << std::endl;
        }
        std::string chrom = header->target_name[record->core.tid];
        int32_t startPos = record->core.pos;  // 0-based start position
        int32_t endPos = bam_endpos(record);  // 1-based end position

        int k = g.getK();
        std::vector<Node *> nodesInRange;
        std::vector<uint> SNVsInRange = extractRange(SNVs, startPos, endPos);
        std::vector<uint> k_SNVs;

        for (int i = 0; i < SNVsInRange.size() - k + 1 && SNVsInRange.size() > k; i++) {
            k_SNVs.assign(SNVsInRange.begin() + i, SNVsInRange.begin() + i + k);

            // Get the two nodes corresponding to this SNV position
            std::vector<Node *> SNVnodes = g.getSNVNodes(chrom, k_SNVs);

            nodesInRange.insert(nodesInRange.begin(), SNVnodes.begin(), SNVnodes.end());
        }

        // Convert BAM record sequence to a string
        std::string readSeq = getReadSequence(record);
        
        std::vector<Node *> nodesToConnect;
        for (auto node : nodesInRange) {
            readXnodes++;
            bool allAlsSupported = true;
            for (char i = 0; i < k; i++) {
                uint snvPosition = std::stoul(node->posArr()[i]);
                std::string allele = node->baseArr()[i];

                allAlsSupported = allAlsSupported && doesReadSupportAllele(record, snvPosition, allele);
            }
            if (allAlsSupported) {
                // std::cout << "Position " << snvPosition << " matches the sequence " << alt << " on this read." << std::endl;
                matchedNodes++;
                nodesToConnect.push_back(node);
            }
            else {
                // std::cout << "Position " << snvPosition << " does not match the sequence " << alt << " on this read." << std::endl;
            }
        }
        g.connectNodeGroup(nodesToConnect);
    }
    std::cout << "Number of reads: " << readNo << std::endl;
    std::cout << "Number of matched nodes: " << matchedNodes << "/" << readXnodes << std::endl;
    // Clean up
    bam_destroy1(record);
    bam_hdr_destroy(header);
    sam_close(inFile);
}


int main(int argc, char* argv[]) {
    std::vector <uint> SNVs;
    Graph g = Graph(1);
    parse_vcf_file(argv[1], g, SNVs);
    connectGraphNodes(argv[2], SNVs, g);

    g.exportToDot("graph.dot");

    return 0;
}