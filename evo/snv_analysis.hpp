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
    float vaf = 0.0f;
    bool needs_read_vaf = false;
    std::string hp_label = "UNKNOWN";
    float hp_proportion = 0.0f;
    int hp1_ref_reads = 0;
    int hp1_alt_reads = 0;
    int hp2_ref_reads = 0;
    int hp2_alt_reads = 0;
    int nohp_alt_reads = 0;
    int total_coverage_reads = 0;
};

struct PositionInfo {
   char ref;
   char alt;
   std::unordered_map<std::string, char> readBase;
   std::unordered_map<std::string, int> readHP; // read name -> HP tag (1/2)
};

struct GermlineSNV {
    std::string chrom;
    int pos;
    char ref;
    char alt;
    int alt_haplotype;
};

struct GermlineVcfStats {
    int total_records = 0;
    int pass_records = 0;
    int eligible_records = 0;
};

struct GermlineVote {
    int hp1 = 0;
    int hp2 = 0;
};

// Parse the somatic VCF and enforce unambiguous sample selection.
bool readVCF(
    const std::string& vcfFile,
    const std::string& requestedSampleName,
    std::vector<SNV>& snvs,
    std::string& resolvedSampleName
);

bool readGermlineVCF(
    const std::string& vcfFile,
    const std::string& sampleName,
    std::vector<GermlineSNV>& snvs,
    GermlineVcfStats& stats
);

// Parse mpileup and attach BAM or externally inferred read haplotypes.
std::unordered_map<std::string, PositionInfo> readMpileup(
    const std::string& mpileupFile,
    const std::vector<SNV>& snvs,
    bool parseBamHpTags = true,
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

// Populate per-SNV HP label/proportion/read counts from pileup info.
void annotateSNVHpSummary(SNV& snv, const PositionInfo& p);

// Calculate ALT / (REF + ALT) from unique pileup reads.
// Returns false when no reads support either allele.
bool calculateVafFromReads(
    const SNV& snv,
    const PositionInfo& p,
    float& vaf
);

// Compare two SNVs and print their co-occurrence counts
void compareSNVs(const SNV& s1, const SNV& s2,
                 const PositionInfo& p1, const PositionInfo& p2,
                 std::ostream& out);

#endif // SNV_ANALYSIS_HPP
