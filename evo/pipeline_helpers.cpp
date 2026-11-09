#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <unordered_map>
#include <vector>
#include <filesystem>
#include <algorithm>
#include <cstdint>
#include <sstream>

#include "pipeline.hpp"

namespace fs = std::filesystem;

static void printUsage(const char* executable) {
    std::cerr << "Internal chromosome worker; use run_genome_by_chr.py for original inputs.\n";
    std::cerr << "Usage: " << executable
              << " --somatic-vcf <VCF_FILE> --bam <BAM_FILE> [--bam <BAM_FILE> ...] --output-dir <OUTPUT_DIR>"
              << " [--window-size <2-5>] [--window-stride <N>]"
              << " [--germline-vcf <VCF_FILE>]"
              << " [--exclude-regions-bed <BED_FILE>]..."
              << " [--divergent-same-hp]\n"
              << "Relationship inference: shared trees (default: 4 mutations, stride 3); "
              << "one-read-per-required-category graph support gate.\n"
              << "--divergent-same-hp requires --germline-vcf.\n";
}

ArgumentResult parsePipelineArguments(int argc, char* argv[], PipelineOptions& options) {
    options.executablePath = argv[0];
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            printUsage(argv[0]);
            return ArgumentResult::Help;
        } else if (arg == "--somatic-vcf" && i + 1 < argc) {
            options.vcfFile = argv[++i];
        } else if (arg == "--bam" && i + 1 < argc) {
            options.bamFiles.push_back(argv[++i]);
        } else if (arg == "--output-dir" && i + 1 < argc) {
            options.outDir = argv[++i];
        } else if (arg == "--window-size" && i + 1 < argc) {
            options.windowSize = std::stoi(argv[++i]);
        } else if (arg == "--window-stride" && i + 1 < argc) {
            options.windowStride = std::stoi(argv[++i]);
        } else if (arg == "--germline-vcf" && i + 1 < argc) {
            options.germlineVcfFile = argv[++i];
        } else if (arg == "--exclude-regions-bed" && i + 1 < argc) {
            options.excludeRegionsBedFiles.push_back(argv[++i]);
        } else if (arg == "--divergent-same-hp") {
            options.divergentSameHp = true;
        } else if (arg.rfind("--", 0) == 0) {
            std::cerr << "Error: unknown or incomplete argument: "
                      << arg << "\n";
            printUsage(argv[0]);
            return ArgumentResult::Error;
        } else {
            std::cerr << "Error: positional argument not supported: "
                      << arg << "\n";
            printUsage(argv[0]);
            return ArgumentResult::Error;
        }
    }

    if (options.vcfFile.empty() || options.bamFiles.empty() || options.outDir.empty()) {
        std::cerr << "Error: missing required arguments.\n";
        printUsage(argv[0]);
        return ArgumentResult::Error;
    }
    if (options.divergentSameHp && options.germlineVcfFile.empty()) {
        std::cerr << "Error: --divergent-same-hp requires --germline-vcf.\n";
        return ArgumentResult::Error;
    }
    if (!fs::exists(options.vcfFile)) {
        std::cerr << "Error: VCF file not found: " << options.vcfFile << "\n";
        return ArgumentResult::Error;
    }
    for (const auto& bamFile : options.bamFiles) {
        if (!fs::exists(bamFile)) {
            std::cerr << "Error: BAM file not found: " << bamFile << "\n";
            return ArgumentResult::Error;
        }
    }
    if (!options.germlineVcfFile.empty() && !fs::exists(options.germlineVcfFile)) {
        std::cerr << "Error: germline VCF file not found: "
                  << options.germlineVcfFile << "\n";
        return ArgumentResult::Error;
    }
    return ArgumentResult::Run;
}


struct BedInterval {
    std::int64_t start;
    std::int64_t end;
};

using BedIntervalsByChrom =
    std::unordered_map<std::string, std::vector<BedInterval>>;


static bool readExclusionBeds(
    const std::vector<std::string>& bedFiles,
    BedIntervalsByChrom& intervalsByChrom,
    std::size_t& intervalCount
) {
    intervalsByChrom.clear();
    intervalCount = 0;

    for (const auto& bedFile : bedFiles) {
        std::ifstream in(bedFile);
        if (!in.is_open()) {
            std::cerr << "Error opening exclusion BED file: " << bedFile << "\n";
            return false;
        }

        std::string line;
        while (std::getline(in, line)) {
            const std::size_t first = line.find_first_not_of(" \t\r\n");
            if (first == std::string::npos || line[first] == '#') continue;

            std::string chrom;
            std::int64_t start = 0;
            std::int64_t end = 0;
            std::istringstream fields(line);
            fields >> chrom >> start >> end;
            intervalsByChrom[chrom].push_back({start, end});
            intervalCount++;
        }
        if (in.bad()) {
            std::cerr << "Error reading exclusion BED: " << bedFile << "\n";
            return false;
        }
    }

    for (auto& [chrom, intervals] : intervalsByChrom) {
        std::sort(
            intervals.begin(),
            intervals.end(),
            [](const BedInterval& a, const BedInterval& b) {
                if (a.start != b.start) return a.start < b.start;
                return a.end < b.end;
            }
        );
        std::vector<BedInterval> merged;
        merged.reserve(intervals.size());
        for (const BedInterval& interval : intervals) {
            if (merged.empty() || interval.start > merged.back().end) {
                merged.push_back(interval);
            } else {
                merged.back().end = std::max(merged.back().end, interval.end);
            }
        }
        intervals = std::move(merged);
    }
    return true;
}


static bool somaticSnvIsExcluded(
    const SNV& snv,
    const BedIntervalsByChrom& intervalsByChrom
) {
    const auto chromIt = intervalsByChrom.find(snv.chrom);
    if (chromIt == intervalsByChrom.end() || chromIt->second.empty()) {
        return false;
    }

    // BED is 0-based and half-open; VCF/SNV positions are 1-based.
    const std::int64_t position0 = static_cast<std::int64_t>(snv.pos) - 1;
    const auto& intervals = chromIt->second;
    auto intervalIt = std::upper_bound(
        intervals.begin(),
        intervals.end(),
        position0,
        [](std::int64_t position, const BedInterval& interval) {
            return position < interval.start;
        }
    );
    if (intervalIt == intervals.begin()) return false;
    --intervalIt;
    return position0 < intervalIt->end;
}

std::string shellQuote(const std::string& value) {
    std::string quoted = "'";
    for (char c : value) {
        quoted += (c == '\'') ? "'\\''" : std::string(1, c);
    }
    return quoted + "'";
}

bool generateGraphs(const PipelineOptions& options, int totalSnvs) {
    const fs::path script = fs::absolute(options.executablePath).parent_path() / "graph_ops/main.py";
    const fs::path chunkDir = fs::path(options.outDir) / "chunk_files";
    const fs::path graphDir = fs::path(options.outDir) / "graphs";
    const std::string command = "python " + shellQuote(script.string()) + " " +
        shellQuote(chunkDir.string()) + " --outdir " + shellQuote(graphDir.string()) +
        " --total-snvs " + std::to_string(totalSnvs);

    std::cout << "\n=== Generating Cytoscape graphs ===\n";
    std::cout << "Running command: " << command << std::endl;
    if (std::system(command.c_str()) != 0) {
        std::cerr << "Error: graph generation failed\n";
        return false;
    }
    std::cout << "Graphs generated in " << graphDir << "\n";
    return true;
}

std::string hpAssignmentFromCounts(
    int hp1,
    int hp2,
    int minHaplotaggedAltReads
) {
    const int total = hp1 + hp2;
    if (total < minHaplotaggedAltReads) return "UNKNOWN";
    const int minor = std::min(hp1, hp2);
    if (100LL * minor > 15LL * total) return "MIXED";
    if (hp1 > hp2) return "HP1";
    if (hp2 > hp1) return "HP2";
    return "MIXED";
}

std::string haplotypeVaf(int ref_reads, int alt_reads) {
    float vaf = 0.0f;
    if (!calculateHaplotypeVaf(ref_reads, alt_reads, vaf)) return "NA";
    return std::to_string(vaf);
}

void writeSnvCountsTsv(
    const fs::path& outDir,
    const std::vector<SNV>& snvs,
    const std::unordered_map<std::string, HaplotypeSnvSummary>& by_key
) {
    const fs::path out_path = outDir / "snv_counts.tsv";
    std::ofstream out(out_path);
    if (!out.is_open()) {
        throw std::runtime_error("Cannot write mutation summaries: " + out_path.string());
    }

    out << "CHR\tPOSITION\tTOTAL_COVERAGE"
        << "\tHP1_SUPPORTING\tHP1_UNSUPPORTING"
        << "\tHP2_SUPPORTING\tHP2_UNSUPPORTING"
        << "\tUNKNOWN_SUPPORTING\tUNKNOWN_UNSUPPORTING"
        << "\tHP1_HVAF\tHP2_HVAF\tHP_ASSIGNMENT\n";
    for (const auto& s : snvs) {
        const std::string key = s.chrom + ":" + std::to_string(s.pos);
        const auto& summary = by_key.at(key);
        out << s.chrom << '\t'
            << s.pos << '\t'
            << summary.total_coverage << '\t'
            << summary.hp1_alt_reads << '\t'
            << summary.hp1_ref_reads << '\t'
            << summary.hp2_alt_reads << '\t'
            << summary.hp2_ref_reads << '\t'
            << summary.nohp_alt_reads << '\t'
            << summary.nohp_ref_reads << '\t'
            << haplotypeVaf(summary.hp1_ref_reads, summary.hp1_alt_reads) << '\t'
            << haplotypeVaf(summary.hp2_ref_reads, summary.hp2_alt_reads) << '\t'
            << hpAssignmentFromCounts(
                   summary.hp1_alt_reads, summary.hp2_alt_reads, 1
               ) << '\n';
    }
    out.close();
    if (!out) throw std::runtime_error("Cannot finish mutation summaries: " + out_path.string());
}


std::string mpileupCommand(
    const std::string& bamArguments, const std::string& region,
    const fs::path& positionsFile, const fs::path& outputFile
) {
    return "samtools mpileup -q 0 -Q 0 --output-QNAME -r " + shellQuote(region) + " -l " + shellQuote(positionsFile.string()) +
           bamArguments + " -o " + shellQuote(outputFile.string());
}


bool filterExcludedRegions(
    std::vector<SNV>& snvs, const std::vector<std::string>& excludeRegionsBedFiles
) {
    if (excludeRegionsBedFiles.empty()) return true;
    BedIntervalsByChrom intervalsByChrom;
    std::size_t bedIntervalCount = 0;
    if (!readExclusionBeds(
            excludeRegionsBedFiles,
            intervalsByChrom,
            bedIntervalCount)) {
        return false;
    }
    const std::size_t inputSomaticSnvs = snvs.size();
    snvs.erase(
        std::remove_if(
            snvs.begin(),
            snvs.end(),
            [&intervalsByChrom](const SNV& snv) {
                return somaticSnvIsExcluded(snv, intervalsByChrom);
            }
        ),
        snvs.end()
    );
    const std::size_t excludedSomaticSnvs = inputSomaticSnvs - snvs.size();
    std::cout << "[somatic_region_filter] bed=";
    for (std::size_t i = 0; i < excludeRegionsBedFiles.size(); ++i) {
        if (i > 0) std::cout << ",";
        std::cout << excludeRegionsBedFiles[i];
    }
    std::cout << " bed_files=" << excludeRegionsBedFiles.size()
              << " bed_intervals=" << bedIntervalCount
              << " input_snvs=" << inputSomaticSnvs
              << " excluded_snvs=" << excludedSomaticSnvs
              << " retained_snvs=" << snvs.size() << "\n";
    return true;
}

bool writeReadHaplotypeVotes(
    const fs::path& readVotePath, std::vector<ReadHaplotypeVoteRow>& readVoteRows
) {
    std::sort(
        readVoteRows.begin(),
        readVoteRows.end(),
        [](const ReadHaplotypeVoteRow& a,
           const ReadHaplotypeVoteRow& b) {
            if (a.chrom != b.chrom) return a.chrom < b.chrom;
            return a.read_name < b.read_name;
        }
    );
    std::ofstream readVoteOut(readVotePath);
    if (!readVoteOut.is_open()) {
        std::cerr << "Error: cannot write " << readVotePath << "\n";
        return false;
    }
    readVoteOut
        << "CHROM\tREAD_NAME\tHP1_GERMLINE_VOTES"
        << "\tHP2_GERMLINE_VOTES"
        << "\tTOTAL_INFORMATIVE_GERMLINE_MUTATIONS"
        << "\tMINORITY_FRACTION\tASSIGNMENT\n";
    for (const auto& row : readVoteRows) {
        readVoteOut << row.chrom << "\t"
                    << row.read_name << "\t"
                    << row.hp1_votes << "\t"
                    << row.hp2_votes << "\t"
                    << (row.hp1_votes + row.hp2_votes) << "\t";
        if (row.minority_fraction < 0.0) {
            readVoteOut << "NA";
        } else {
            readVoteOut << row.minority_fraction;
        }
        readVoteOut << "\t" << row.assignment << "\n";
    }
    readVoteOut.close();
    if (!readVoteOut) {
        std::cerr << "Error: cannot finish " << readVotePath << "\n";
        return false;
    }
    return true;
}
