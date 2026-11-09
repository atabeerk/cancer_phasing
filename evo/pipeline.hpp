#ifndef PIPELINE_HPP
#define PIPELINE_HPP

#include <filesystem>
#include <fstream>
#include <string>
#include <stdexcept>
#include <vector>

#include "snv_analysis.hpp"

struct PipelineOptions {
    std::string executablePath;
    std::string vcfFile;
    std::vector<std::string> bamFiles;
    std::string outDir;
    std::string germlineVcfFile;
    std::vector<std::string> excludeRegionsBedFiles;
    int windowSize = 4;
    int windowStride = 3;
    bool divergentSameHp = false;
};

enum class ArgumentResult { Run, Help, Error };

ArgumentResult parsePipelineArguments(int argc, char* argv[], PipelineOptions& options);
int runChromosomePipeline(const PipelineOptions& options);

struct ReadHaplotypeVoteRow {
    std::string chrom;
    std::string read_name;
    int hp1_votes;
    int hp2_votes;
    double minority_fraction;
    std::string assignment;
};

bool filterExcludedRegions(
    std::vector<SNV>& snvs, const std::vector<std::string>& bedFiles
);
std::string hpAssignmentFromCounts(int hp1, int hp2, int minHaplotaggedAltReads);
void writeSnvCountsTsv(
    const std::filesystem::path& outDir,
    const std::vector<SNV>& snvs,
    const std::unordered_map<std::string, HaplotypeSnvSummary>& summaries
);
bool writeReadHaplotypeVotes(
    const std::filesystem::path& path, std::vector<ReadHaplotypeVoteRow>& rows
);
std::string mpileupCommand(
    const std::string& bamArguments, const std::string& region,
    const std::filesystem::path& positionsFile, const std::filesystem::path& outputFile
);
std::string shellQuote(const std::string& value);
bool generateGraphs(const PipelineOptions& options, int totalSnvs);

// Both somatic and germline pileups use 0-based, half-open position intervals.
template <typename Variant>
void writePileupPositions(
    const std::filesystem::path& path, const std::vector<Variant>& variants
) {
    std::ofstream out(path);
    for (const auto& variant : variants) {
        out << variant.chrom << "\t" << (variant.pos - 1)
            << "\t" << variant.pos << "\n";
    }
    out.close();
    if (!out) throw std::runtime_error("Cannot write pileup positions: " + path.string());
}

#endif
