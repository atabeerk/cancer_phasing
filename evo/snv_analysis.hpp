#ifndef SNV_ANALYSIS_HPP
#define SNV_ANALYSIS_HPP

#include <string>
#include <vector>
#include <unordered_map>

// --- Structs ---

struct SNV {
    std::string chrom;
    int pos;
    char ref;
    char alt;
    float vaf;
};

struct PositionInfo {
   char ref;
   char alt;
   std::unordered_map<std::string, char> readBase;
};

// Parse VCF (supports .vcf and .vcf.gz via htslib)
std::vector<SNV> readVCF(const std::string& vcfFile);

// Parse mpileup (-s) file (7 columns: chrom pos ref depth bases qualities readnames)
std::unordered_map<std::string, PositionInfo> readMpileup(const std::string& mpileupFile, const std::vector<SNV>& snvs);

// Compare two SNVs and print their co-occurrence counts
void compareSNVs(const SNV& s1, const SNV& s2,
                 const PositionInfo& p1, const PositionInfo& p2,
                 std::ostream& out);

#endif // SNV_ANALYSIS_HPP
