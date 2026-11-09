#ifndef SNV_ANALYSIS_HPP
#define SNV_ANALYSIS_HPP

#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>

// --- Structs ---

struct SNV {
    std::string chrom;
    int pos;
    char ref;
    char alt;
};

struct HaplotypeSnvSummary {
    int total_coverage = 0;
    int hp1_ref_reads = 0;
    int hp1_alt_reads = 0;
    int hp2_ref_reads = 0;
    int hp2_alt_reads = 0;
    int nohp_ref_reads = 0;
    int nohp_alt_reads = 0;

    int altReadCount() const {
        return hp1_alt_reads + hp2_alt_reads + nohp_alt_reads;
    }
    int alleleDepth() const {
        return altReadCount() + hp1_ref_reads + hp2_ref_reads + nohp_ref_reads;
    }
};

struct PositionInfo {
   char ref;
   char alt;
   std::unordered_map<std::string, char> readBase;
};

struct GermlineSNV {
    std::string chrom;
    int pos;
    char ref;
    char alt;
    int alt_haplotype;
};

struct GermlineVote {
    int hp1 = 0;
    int hp2 = 0;
};

// Read chromosome VCFs prepared and filtered by Python preprocessing.
bool readVCF(const std::string& vcfFile, std::vector<SNV>& snvs);
bool readGermlineVCF(const std::string& vcfFile, std::vector<GermlineSNV>& snvs);

// Stream positions and count alleles using optional germline-derived haplotypes.
std::unordered_map<std::string, HaplotypeSnvSummary> summarizeSomaticMpileup(
    const std::string& mpileupFile,
    const std::vector<SNV>& snvs,
    const std::unordered_map<std::string, int>* externalReadHp = nullptr
);

void collectReadNamesFromMpileup(
    const std::string& mpileupFile,
    std::unordered_set<std::string>& readNames
);

std::size_t accumulateGermlineVotes(
    const std::string& mpileupFile,
    const std::vector<GermlineSNV>& snvs,
    const std::unordered_set<std::string>& candidateReads,
    std::unordered_map<std::string, GermlineVote>& votes
);

// Count unique REF/ALT reads by haplotype for one somatic mutation.
HaplotypeSnvSummary summarizeSnvReads(
    const PositionInfo& p, const std::unordered_map<std::string, int>* readHaplotypes
);

// Calculate ALT / (REF + ALT) for one haplotype.
// Returns false when that haplotype has no REF or ALT reads.
bool calculateHaplotypeVaf(
    int refReads,
    int altReads,
    float& vaf
);

#endif // SNV_ANALYSIS_HPP
