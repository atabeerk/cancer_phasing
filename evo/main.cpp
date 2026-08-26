#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <filesystem>
#include <algorithm>
#include <chrono>
#include <random>
#include <cstdint>
#include <sstream>

#include "htslib/vcf.h"
#include "htslib/sam.h"

#include "snv_analysis.hpp"
#include "process_snv_output.hpp"

namespace fs = std::filesystem;


bool commandAvailable(const std::string& cmd) {
    const std::string probe = "command -v " + cmd + " >/dev/null 2>&1";
    return std::system(probe.c_str()) == 0;
}


struct BedInterval {
    std::int64_t start;
    std::int64_t end;
};

using BedIntervalsByChrom =
    std::unordered_map<std::string, std::vector<BedInterval>>;


bool readExclusionBed(
    const std::string& bedFile,
    BedIntervalsByChrom& intervalsByChrom,
    std::size_t& intervalCount
) {
    intervalsByChrom.clear();
    intervalCount = 0;

    std::ifstream in(bedFile);
    if (!in.is_open()) {
        std::cerr << "Error opening exclusion BED file: " << bedFile << "\n";
        return false;
    }

    std::string line;
    std::size_t lineNumber = 0;
    while (std::getline(in, line)) {
        lineNumber++;
        const std::size_t first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos || line[first] == '#') continue;

        std::string chrom;
        std::int64_t start = 0;
        std::int64_t end = 0;
        std::istringstream fields(line);
        if (!(fields >> chrom >> start >> end) || start < 0 || end <= start) {
            std::cerr << "Error: invalid BED interval at " << bedFile
                      << ":" << lineNumber << "\n";
            return false;
        }
        intervalsByChrom[chrom].push_back({start, end});
        intervalCount++;
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


bool somaticSnvIsExcluded(
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

int detectMaxPairDistanceFromBams(const std::vector<std::string>& bamFiles) {
    const auto t0 = std::chrono::steady_clock::now();
    int max_ref_span = 0;
    for (const auto& bamFile : bamFiles) {
    samFile* fp = sam_open(bamFile.c_str(), "r");
    if (!fp) {
        const auto t1 = std::chrono::steady_clock::now();
        const double elapsed_s =
            std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0).count();
        std::cerr << "Warning: cannot open BAM for read-length scan: " << bamFile
                  << ". Falling back to 200000 bp max pair distance.\n";
        std::cout << "[max_span_scan] elapsed_seconds=" << elapsed_s
                  << " fallback_max_pair_distance_bp=200000\n";
        return 200000;
    }

    bam_hdr_t* hdr = sam_hdr_read(fp);
    if (!hdr) {
        const auto t1 = std::chrono::steady_clock::now();
        const double elapsed_s =
            std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0).count();
        std::cerr << "Warning: cannot read BAM header for read-length scan: " << bamFile
                  << ". Falling back to 200000 bp max pair distance.\n";
        sam_close(fp);
        std::cout << "[max_span_scan] elapsed_seconds=" << elapsed_s
                  << " fallback_max_pair_distance_bp=200000\n";
        return 200000;
    }

    bam1_t* rec = bam_init1();
    if (!rec) {
        const auto t1 = std::chrono::steady_clock::now();
        const double elapsed_s =
            std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0).count();
        std::cerr << "Warning: cannot allocate BAM record for read-length scan. "
                  << "Falling back to 200000 bp max pair distance.\n";
        bam_hdr_destroy(hdr);
        sam_close(fp);
        std::cout << "[max_span_scan] elapsed_seconds=" << elapsed_s
                  << " fallback_max_pair_distance_bp=200000\n";
        return 200000;
    }

    while (sam_read1(fp, hdr, rec) >= 0) {
        if (rec->core.flag & BAM_FUNMAP) {
            continue;
        }
        if (rec->core.n_cigar <= 0) {
            continue;
        }
        int ref_span = bam_cigar2rlen(rec->core.n_cigar, bam_get_cigar(rec));
        max_ref_span = std::max(max_ref_span, ref_span);
    }

    bam_destroy1(rec);
    bam_hdr_destroy(hdr);
    sam_close(fp);

    }
    if (max_ref_span <= 0) {
        const auto t1 = std::chrono::steady_clock::now();
        const double elapsed_s =
            std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0).count();
        std::cerr << "Warning: BAM aligned-span scan found no usable mapped alignments. "
                  << "Falling back to 200000 bp max pair distance.\n";
        std::cout << "[max_span_scan] elapsed_seconds=" << elapsed_s
                  << " fallback_max_pair_distance_bp=200000\n";
        return 200000;
    }

    // Round up to the next kb, e.g. 200341 -> 201000.
    int rounded_kb = ((max_ref_span + 999) / 1000) * 1000;
    const auto t1 = std::chrono::steady_clock::now();
    const double elapsed_s =
        std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0).count();
    std::cout << "[max_span_scan] elapsed_seconds=" << elapsed_s
              << " max_aligned_ref_span_bp=" << max_ref_span
              << " selected_max_pair_distance_bp=" << rounded_kb << "\n";
    return rounded_kb;
}

bool detectBamsHaplotagged(
    const std::vector<std::string>& bamFiles,
    int sampleReads = 1000
) {
    const auto t0 = std::chrono::steady_clock::now();
    std::vector<uint8_t> sampled_has_hp;
    sampled_has_hp.reserve(sampleReads > 0 ? sampleReads : 0);
    std::uint64_t total_seen = 0;
    std::mt19937_64 rng(std::random_device{}());
    for (const auto& bamFile : bamFiles) {
    samFile* fp = sam_open(bamFile.c_str(), "r");
    if (!fp) {
        std::cerr << "Warning: cannot open BAM for haplotag detection: " << bamFile << "\n";
        return false;
    }
    bam_hdr_t* hdr = sam_hdr_read(fp);
    if (!hdr) {
        std::cerr << "Warning: cannot read BAM header for haplotag detection: " << bamFile << "\n";
        sam_close(fp);
        return false;
    }
    bam1_t* rec = bam_init1();
    if (!rec) {
        std::cerr << "Warning: cannot allocate BAM record for haplotag detection.\n";
        bam_hdr_destroy(hdr);
        sam_close(fp);
        return false;
    }


    while (sam_read1(fp, hdr, rec) >= 0) {
        ++total_seen;
        const bool has_hp = (bam_aux_get(rec, "HP") != nullptr);
        if (static_cast<int>(sampled_has_hp.size()) < sampleReads) {
            sampled_has_hp.push_back(has_hp ? 1 : 0);
        } else if (sampleReads > 0) {
            std::uniform_int_distribution<std::uint64_t> dist(0, total_seen - 1);
            const std::uint64_t j = dist(rng);
            if (j < static_cast<std::uint64_t>(sampleReads)) {
                sampled_has_hp[static_cast<std::size_t>(j)] = has_hp ? 1 : 0;
            }
        }
    }

    bam_destroy1(rec);
    bam_hdr_destroy(hdr);
    sam_close(fp);
    }
    int hp_tagged_sampled_reads = 0;
    for (uint8_t v : sampled_has_hp) {
        if (v != 0) hp_tagged_sampled_reads++;
    }
    const bool has_hp_tag = (hp_tagged_sampled_reads > 0);


    std::cout << "[haplotag_detect] sampled_reads=" << sampled_has_hp.size()
              << " total_seen=" << total_seen
              << " hp_tagged_sampled_reads=" << hp_tagged_sampled_reads
              << " has_hp_tag=" << (has_hp_tag ? 1 : 0) << "\n";
    const double elapsed_s =
        std::chrono::duration_cast<std::chrono::duration<double>>(
            std::chrono::steady_clock::now() - t0
        ).count();
    std::cout << "[haplotag_timing] mode=bam_hp_tags"
              << " stage=bam_hp_detection seconds="
              << elapsed_s << "\n";
    return has_hp_tag;
}

void vcfToBed(const std::string& vcfFile, const std::string& bedFile) {
    htsFile* fp = bcf_open(vcfFile.c_str(), "r");
    if (!fp) {
        std::cerr << "Error: cannot open VCF file: " << vcfFile << "\n";
        return;
    }

    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    if (!hdr) {
        std::cerr << "Error: cannot read VCF header\n";
        bcf_close(fp);
        return;
    }

    bcf1_t* rec = bcf_init();
    if (!rec) {
        std::cerr << "Error: cannot allocate bcf1_t\n";
        bcf_hdr_destroy(hdr);
        bcf_close(fp);
        return;
    }

    std::ofstream bed(bedFile);
    if (!bed.is_open()) {
        std::cerr << "Error: cannot open BED file for writing: " << bedFile << "\n";
        bcf_destroy(rec);
        bcf_hdr_destroy(hdr);
        bcf_close(fp);
        return;
    }

    while (bcf_read(fp, hdr, rec) == 0) {
        bcf_unpack(rec, BCF_UN_STR);
        const char* chrom = bcf_hdr_id2name(hdr, rec->rid);
        int start = rec->pos;
        int end = rec->pos + 1;
        bed << chrom << "\t" << start << "\t" << end << "\n";
    }

    bed.close();
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(fp);

    std::cout << "BED file written to: " << bedFile << "\n";
}


void generateGraphs(const std::string& outDir, int totalSnvs) {
    std::string chunkDir = fs::path(outDir) / "chunk_files";
    std::string graphDir = fs::path(outDir) / "graphs";
    std::string cmd = "python graph_ops/main.py " + chunkDir + " --outdir " + graphDir
                    + " --total-snvs " + std::to_string(totalSnvs);

    std::cout << "\n=== Generating Cytoscape graphs ===\n";
    std::cout << "Running command: " << cmd << "\n";

    int ret = system(cmd.c_str());
    if (ret != 0) {
        std::cerr << "Warning: Graph generation failed with exit code " << ret << "\n";
    } else {
        std::cout << "Graphs generated in " << graphDir << "\n";
    }
}

struct HaplotypeSnvSummary {
    int total_coverage = 0;
    int hp1_ref_reads = 0;
    int hp1_alt_reads = 0;
    int hp2_ref_reads = 0;
    int hp2_alt_reads = 0;
    int nohp_alt_reads = 0;
};

static std::string hpAssignmentFromCounts(
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

static std::string haplotypeVaf(int ref_reads, int alt_reads) {
    const int allele_depth = ref_reads + alt_reads;
    if (allele_depth == 0) return "NA";
    return std::to_string(static_cast<double>(alt_reads) /
                          static_cast<double>(allele_depth));
}

static void writeHaplotaggedSnvsTsv(
    const fs::path& outDir,
    const std::vector<SNV>& snvs,
    const std::unordered_map<std::string, HaplotypeSnvSummary>& by_key,
    const std::unordered_set<std::string>& skippedSnvKeys,
    int minHaplotaggedAltReads
) {
    const fs::path out_path = outDir / "haplotagged_snvs.tsv";
    std::ofstream out(out_path);
    if (!out.is_open()) {
        std::cerr << "Warning: failed to write " << out_path << "\n";
        return;
    }

    out << "CHR\tPOSITION\tTOTAL_COVERAGE\tHP_COUNTS"
        << "\tHP1_HVAF\tHP2_HVAF\tHP_ASSIGNMENT\n";
    for (const auto& s : snvs) {
        const std::string key = s.chrom + ":" + std::to_string(s.pos);
        if (skippedSnvKeys.count(key) != 0) continue;
        int total_coverage = 0;
        int hp1_ref = 0;
        int hp1 = 0;
        int hp2_ref = 0;
        int hp2 = 0;
        int nohp = 0;
        auto it = by_key.find(key);
        if (it != by_key.end()) {
            total_coverage = it->second.total_coverage;
            hp1_ref = it->second.hp1_ref_reads;
            hp1 = it->second.hp1_alt_reads;
            hp2_ref = it->second.hp2_ref_reads;
            hp2 = it->second.hp2_alt_reads;
            nohp = it->second.nohp_alt_reads;
        }
        out << s.chrom << '\t'
            << s.pos << '\t'
            << total_coverage << '\t'
            << hp1 << "/" << hp2 << "/" << nohp << '\t'
            << haplotypeVaf(hp1_ref, hp1) << '\t'
            << haplotypeVaf(hp2_ref, hp2) << '\t'
            << hpAssignmentFromCounts(
                   hp1, hp2, minHaplotaggedAltReads
               ) << '\n';
    }
    out.close();
}


struct ChunkWork {
    std::string chrom;
    int region_start;
    int region_end;
    std::vector<SNV> halo_snvs;
    fs::path mpileup_out;
    fs::path chunk_out;
};

struct ReadHaplotypeVoteRow {
    std::string chrom;
    std::string read_name;
    int hp1_votes;
    int hp2_votes;
    double minority_fraction;
    std::string assignment;
};

static constexpr int kMinInformativeGermlineMutations = 2;
static constexpr int kMinBamHpTaggedAltReadsPerSomaticMutation = 2;

static double elapsedSeconds(
    const std::chrono::steady_clock::time_point& start
) {
    return std::chrono::duration_cast<std::chrono::duration<double>>(
        std::chrono::steady_clock::now() - start
    ).count();
}

int main(int argc, char* argv[]) {
    const auto main_t0 = std::chrono::steady_clock::now();
    auto printUsage = [&]() {
        std::cerr << "Usage: " << argv[0]
                  << " --somatic-vcf <VCF_FILE> --bam <BAM_FILE> [--bam <BAM_FILE> ...] --output-dir <OUTPUT_DIR>"
                  << " [--min-reads <N>] [--somatic-sample-name <SAMPLE_NAME>]"
                  << " [--germline-vcf <VCF_FILE>]"
                  << " [--germline-vcf-sample-name <SAMPLE_NAME>]"
                  << " [--exclude-regions-bed <BED_FILE>]"
                  << " [--max-pair-distance-kb <N>]"
                  << " [--divergent-same-hp]"
                  << " [--bam-haplotagged <0|1>]\n";
    };

    std::string vcfFile;
    std::vector<std::string> bamFiles;
    std::string outDir;
    std::string germlineVcfFile;
    std::string germlineVcfSampleName;
    std::string excludeRegionsBedFile;
    int minReads = 2;
    std::string vcfSampleName;
    int maxPairDistanceKbArg = -1;
    int bamHaplotaggedArg = -1;
    bool divergentSameHp = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            printUsage();
            return 0;
        } else if (arg == "--somatic-vcf" && i + 1 < argc) {
            vcfFile = argv[++i];
        } else if (arg == "--bam" && i + 1 < argc) {
            bamFiles.push_back(argv[++i]);
        } else if (arg == "--output-dir" && i + 1 < argc) {
            outDir = argv[++i];
        } else if (arg == "--min-reads" && i + 1 < argc) {
            minReads = std::stoi(argv[++i]);
        } else if (arg == "--somatic-sample-name" && i + 1 < argc) {
            vcfSampleName = argv[++i];
        } else if (arg == "--germline-vcf" && i + 1 < argc) {
            germlineVcfFile = argv[++i];
        } else if (arg == "--germline-vcf-sample-name" && i + 1 < argc) {
            germlineVcfSampleName = argv[++i];
        } else if (arg == "--exclude-regions-bed" && i + 1 < argc) {
            excludeRegionsBedFile = argv[++i];
        } else if (arg == "--max-pair-distance-kb" && i + 1 < argc) {
            maxPairDistanceKbArg = std::stoi(argv[++i]);
            if (maxPairDistanceKbArg <= 0) {
                std::cerr << "Error: --max-pair-distance-kb must be > 0\n";
                return 1;
            }
        } else if (arg == "--divergent-same-hp") {
            divergentSameHp = true;
        } else if (arg == "--bam-haplotagged" && i + 1 < argc) {
            bamHaplotaggedArg = std::stoi(argv[++i]);
            if (bamHaplotaggedArg != 0 && bamHaplotaggedArg != 1) {
                std::cerr << "Error: --bam-haplotagged must be 0 or 1\n";
                return 1;
            }
        } else if (arg.rfind("--", 0) == 0) {
            std::cerr << "Error: unknown or incomplete argument: "
                      << arg << "\n";
            printUsage();
            return 1;
        } else {
            std::cerr << "Error: positional argument not supported: "
                      << arg << "\n";
            printUsage();
            return 1;
        }
    }

    if (vcfFile.empty() || bamFiles.empty() || outDir.empty()) {
        std::cerr << "Error: missing required arguments.\n";
        printUsage();
        return 1;
    }
    if (!germlineVcfFile.empty() && germlineVcfSampleName.empty()) {
        std::cerr << "Error: --germline-vcf-sample-name is required with "
                  << "--germline-vcf.\n";
        return 1;
    }
    if (germlineVcfFile.empty() && !germlineVcfSampleName.empty()) {
        std::cerr << "Error: --germline-vcf-sample-name requires "
                  << "--germline-vcf.\n";
        return 1;
    }

    if (!commandAvailable("samtools")) {
        std::cerr << "Error: required executable not found on PATH: samtools\n";
        return 1;
    }
    if (!fs::exists(vcfFile)) {
        std::cerr << "Error: VCF file not found: " << vcfFile << "\n";
        return 1;
    }
    for (const auto& bamFile : bamFiles) {
        if (!fs::exists(bamFile)) {
            std::cerr << "Error: BAM file not found: " << bamFile << "\n";
            return 1;
        }
    }
    std::string bamArguments;
    for (const auto& bamFile : bamFiles) {
        bamArguments += " " + bamFile;
    }
    if (!germlineVcfFile.empty() && !fs::exists(germlineVcfFile)) {
        std::cerr << "Error: germline VCF file not found: "
                  << germlineVcfFile << "\n";
        return 1;
    }
    if (!excludeRegionsBedFile.empty() &&
        !fs::exists(excludeRegionsBedFile)) {
        std::cerr << "Error: exclusion BED file not found: "
                  << excludeRegionsBedFile << "\n";
        return 1;
    }

    const bool useGermlineVcf = !germlineVcfFile.empty();
    const int minHaplotaggedAltReadsPerSomaticMutation =
        useGermlineVcf
            ? 1
            : kMinBamHpTaggedAltReadsPerSomaticMutation;
    if (useGermlineVcf && bamHaplotaggedArg != -1) {
        std::cerr << "Warning: --bam-haplotagged is ignored when "
                  << "--germline-vcf is provided.\n";
    }
    std::cout << "Using MIN_READS = " << minReads << "\n";
    std::cout << "[haplotag_mode] mode="
              << (useGermlineVcf ? "germline_vcf" : "bam_hp_tags")
              << " min_haplotagged_alt_reads_per_somatic_mutation="
              << minHaplotaggedAltReadsPerSomaticMutation
              << "\n";
    std::cout << "[relationship_filter] divergent_same_hp="
              << (divergentSameHp ? 1 : 0) << "\n";

    if (!fs::exists(outDir)) {
        fs::create_directories(outDir);
    }
    const fs::path pileupDir = fs::path(outDir) / "pileup_files";
    const fs::path chunkDir = fs::path(outDir) / "chunk_files";
    fs::create_directories(pileupDir);
    fs::create_directories(chunkDir);

    std::vector<SNV> snvs;
    std::string resolvedSomaticSampleName;
    if (!readVCF(
            vcfFile,
            vcfSampleName,
            snvs,
            resolvedSomaticSampleName)) {
        return 1;
    }

    if (!excludeRegionsBedFile.empty()) {
        BedIntervalsByChrom intervalsByChrom;
        std::size_t bedIntervalCount = 0;
        if (!readExclusionBed(
                excludeRegionsBedFile,
                intervalsByChrom,
                bedIntervalCount)) {
            return 1;
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
        std::cout << "[somatic_region_filter] bed="
                  << excludeRegionsBedFile
                  << " bed_intervals=" << bedIntervalCount
                  << " input_snvs=" << inputSomaticSnvs
                  << " excluded_snvs=" << excludedSomaticSnvs
                  << " retained_snvs=" << snvs.size() << "\n";
    }

    std::vector<GermlineSNV> germlineSnvs;
    GermlineVcfStats germlineVcfStats;
    if (useGermlineVcf &&
        !readGermlineVCF(
            germlineVcfFile,
            germlineVcfSampleName,
            germlineSnvs,
            germlineVcfStats)) {
        return 1;
    }

    std::unordered_map<std::string, HaplotypeSnvSummary> hapSummaries;
    hapSummaries.reserve(snvs.size());
    std::unordered_set<std::string> readVafCandidateKeys;
    for (const auto& snv : snvs) {
        const std::string key =
            snv.chrom + ":" + std::to_string(snv.pos);
        hapSummaries.emplace(
            key,
            HaplotypeSnvSummary{}
        );
        if (snv.needs_read_vaf) {
            readVafCandidateKeys.insert(key);
        }
    }

    std::unordered_map<std::string, std::vector<SNV>> chromSnvs;
    for (const auto& snv : snvs) {
        chromSnvs[snv.chrom].push_back(snv);
    }

    int maxPairDistanceBp = 0;
    if (maxPairDistanceKbArg > 0) {
        maxPairDistanceBp = maxPairDistanceKbArg * 1000;
        std::cout << "Using max SNV pair distance = "
                  << maxPairDistanceBp
                  << " bp (from --max-pair-distance-kb)\n";
    } else {
        maxPairDistanceBp = detectMaxPairDistanceFromBams(bamFiles);
        std::cout << "Using max SNV pair distance = "
                  << maxPairDistanceBp
                  << " bp (rounded from max aligned reference span in BAM)\n";
    }

    bool bamIsHaplotagged = false;
    if (!useGermlineVcf) {
        if (bamHaplotaggedArg == 0 || bamHaplotaggedArg == 1) {
            bamIsHaplotagged = bamHaplotaggedArg == 1;
            std::cout << "Using BAM haplotagged status = "
                      << (bamIsHaplotagged ? 1 : 0)
                      << " (from --bam-haplotagged)\n";
        } else {
            bamIsHaplotagged = detectBamsHaplotagged(bamFiles);
            std::cout << "Using BAM haplotagged status = "
                      << (bamIsHaplotagged ? 1 : 0)
                      << " (auto-detected)\n";
        }
        std::cout << "[haplotag_status] mode=bam_hp_tags haplotagged="
                  << (bamIsHaplotagged ? 1 : 0)
                  << " source="
                  << (bamHaplotaggedArg == -1 ? "auto-detected" : "cli")
                  << "\n";
    }

    double somaticPileupSeconds = 0.0;
    double candidateCollectionSeconds = 0.0;
    double germlinePileupSeconds = 0.0;
    double germlineVoteSeconds = 0.0;
    double somaticParseAttachSeconds = 0.0;
    double somaticHpSummarySeconds = 0.0;
    double pairComparisonSeconds = 0.0;

    std::unordered_map<std::string, float> readVafByKey;
    std::unordered_set<std::string> skippedReadVafKeys;

    const int chunkSize = 10'000'000;
    std::vector<ChunkWork> chunks;
    std::unordered_map<
        std::string,
        std::unordered_set<std::string>
    > candidateReadsByChrom;

    // First pass: create each somatic pileup. In germline mode, also collect
    // the read names that cover at least one somatic mutation.
    for (const auto& chromEntry : chromSnvs) {
        const std::string& chrom = chromEntry.first;
        const std::vector<SNV>& snvList = chromEntry.second;
        if (snvList.empty()) continue;

        std::vector<SNV> sortedSnvList = snvList;
        std::sort(
            sortedSnvList.begin(),
            sortedSnvList.end(),
            [](const SNV& a, const SNV& b) { return a.pos < b.pos; }
        );

        const int minPos = sortedSnvList.front().pos;
        const int maxPos = sortedSnvList.back().pos;
        for (int regionStart = minPos;
             regionStart <= maxPos;
             regionStart += chunkSize) {
            const int regionEnd = regionStart + chunkSize - 1;
            const int haloStart =
                std::max(1, regionStart - maxPairDistanceBp);
            const int haloEnd = regionEnd + maxPairDistanceBp;

            std::vector<SNV> haloSnvs;
            std::size_t coreSnvCount = 0;
            for (const auto& snv : sortedSnvList) {
                if (snv.pos >= haloStart && snv.pos <= haloEnd) {
                    haloSnvs.push_back(snv);
                    if (snv.pos >= regionStart && snv.pos <= regionEnd) {
                        coreSnvCount++;
                    }
                }
            }
            if (coreSnvCount == 0) continue;

            std::cout << "[chunk_halo] chrom=" << chrom
                      << " core_start=" << regionStart
                      << " core_end=" << regionEnd
                      << " halo_start=" << haloStart
                      << " halo_end=" << haloEnd
                      << " core_snvs=" << coreSnvCount
                      << " halo_snvs=" << haloSnvs.size() << "\n";

            const fs::path positionsFile =
                pileupDir /
                ("positions_" + chrom + "_" +
                 std::to_string(regionStart) + ".txt");
            {
                std::ofstream positions(positionsFile);
                for (const auto& snv : haloSnvs) {
                    positions << snv.chrom << "\t"
                              << (snv.pos - 1) << "\t"
                              << snv.pos << "\n";
                }
            }

            const fs::path mpileupOut =
                pileupDir /
                ("mpileup_" + chrom + "_" +
                 std::to_string(regionStart) + ".out");
            std::string command =
                "samtools mpileup -q 0 -Q 0 --output-QNAME";
            if (!useGermlineVcf && bamIsHaplotagged) {
                command += " --output-extra HP";
            }
            command += " -r " + chrom + ":" +
                       std::to_string(haloStart) + "-" +
                       std::to_string(haloEnd);
            command += " -l " + positionsFile.string() +
                       bamArguments +
                       " -o " + mpileupOut.string();

            std::cout << "Running command:\n" << command << "\n";
            const auto pileupT0 = std::chrono::steady_clock::now();
            const int pileupStatus = system(command.c_str());
            somaticPileupSeconds += elapsedSeconds(pileupT0);
            if (pileupStatus != 0) {
                std::cerr << "Error: mpileup failed for region "
                          << chrom << ":" << regionStart
                          << "-" << regionEnd << "\n";
                return 1;
            }

            if (useGermlineVcf) {
                const auto candidateT0 =
                    std::chrono::steady_clock::now();
                collectReadNamesFromMpileup(
                    mpileupOut.string(),
                    candidateReadsByChrom[chrom]
                );
                candidateCollectionSeconds +=
                    elapsedSeconds(candidateT0);
            }

            const fs::path chunkOut =
                chunkDir /
                ("chunk_" + chrom + "_" +
                 std::to_string(regionStart) + "_" +
                 std::to_string(regionEnd) + ".txt");
            chunks.push_back(ChunkWork{
                chrom,
                regionStart,
                regionEnd,
                std::move(haloSnvs),
                mpileupOut,
                chunkOut
            });
        }
    }

    std::unordered_map<
        std::string,
        std::unordered_map<std::string, int>
    > inferredReadHpByChrom;
    std::size_t candidateReadCount = 0;
    std::size_t informativeGermlineObservations = 0;
    std::size_t inferredHp1Reads = 0;
    std::size_t inferredHp2Reads = 0;
    std::size_t mixedReads = 0;
    std::vector<ReadHaplotypeVoteRow> readVoteRows;

    if (useGermlineVcf) {
        std::unordered_map<
            std::string,
            std::vector<GermlineSNV>
        > germlineSnvsByChrom;
        for (const auto& snv : germlineSnvs) {
            germlineSnvsByChrom[snv.chrom].push_back(snv);
        }

        for (const auto& candidateEntry : candidateReadsByChrom) {
            const std::string& chrom = candidateEntry.first;
            const auto& candidateReads = candidateEntry.second;
            candidateReadCount += candidateReads.size();

            auto markerIt = germlineSnvsByChrom.find(chrom);
            if (candidateReads.empty()) {
                continue;
            }
            if (markerIt == germlineSnvsByChrom.end() ||
                markerIt->second.empty()) {
                for (const auto& readName : candidateReads) {
                    readVoteRows.push_back({
                        chrom, readName, 0, 0, -1.0, "UNASSIGNED"
                    });
                }
                continue;
            }

            std::vector<GermlineSNV>& markers = markerIt->second;
            std::sort(
                markers.begin(),
                markers.end(),
                [](const GermlineSNV& a, const GermlineSNV& b) {
                    return a.pos < b.pos;
                }
            );

            const fs::path positionsFile =
                pileupDir / ("germline_positions_" + chrom + ".txt");
            {
                std::ofstream positions(positionsFile);
                for (const auto& marker : markers) {
                    positions << marker.chrom << "\t"
                              << (marker.pos - 1) << "\t"
                              << marker.pos << "\n";
                }
            }

            const fs::path mpileupOut =
                pileupDir / ("germline_mpileup_" + chrom + ".out");
            std::string command =
                "samtools mpileup -q 0 -Q 0 --output-QNAME";
            command += " -r " + chrom + ":" +
                       std::to_string(markers.front().pos) + "-" +
                       std::to_string(markers.back().pos);
            command += " -l " + positionsFile.string() +
                       bamArguments +
                       " -o " + mpileupOut.string();

            std::cout << "Running germline command:\n"
                      << command << "\n";
            const auto pileupT0 = std::chrono::steady_clock::now();
            const int pileupStatus = system(command.c_str());
            germlinePileupSeconds += elapsedSeconds(pileupT0);
            if (pileupStatus != 0) {
                std::cerr << "Error: germline mpileup failed for "
                          << chrom << "\n";
                return 1;
            }

            const auto voteT0 = std::chrono::steady_clock::now();
            std::unordered_map<std::string, GermlineVote> votes;
            informativeGermlineObservations +=
                accumulateGermlineVotes(
                    mpileupOut.string(),
                    markers,
                    candidateReads,
                    votes
                );

            auto& inferredReadHp = inferredReadHpByChrom[chrom];
            inferredReadHp.reserve(votes.size());
            for (const auto& readName : candidateReads) {
                int hp1 = 0;
                int hp2 = 0;
                auto voteIt = votes.find(readName);
                if (voteIt != votes.end()) {
                    hp1 = voteIt->second.hp1;
                    hp2 = voteIt->second.hp2;
                }
                const int total = hp1 + hp2;
                if (total < kMinInformativeGermlineMutations) {
                    readVoteRows.push_back({
                        chrom, readName, hp1, hp2, -1.0, "UNASSIGNED"
                    });
                    continue;
                }

                const int minor = std::min(hp1, hp2);
                const double minorityFraction =
                    static_cast<double>(minor) /
                    static_cast<double>(total);
                if (100LL * minor > 15LL * total || hp1 == hp2) {
                    mixedReads++;
                    readVoteRows.push_back({
                        chrom, readName, hp1, hp2,
                        minorityFraction, "MIXED"
                    });
                } else if (hp1 > hp2) {
                    inferredReadHp[readName] = 1;
                    inferredHp1Reads++;
                    readVoteRows.push_back({
                        chrom, readName, hp1, hp2,
                        minorityFraction, "HP1"
                    });
                } else {
                    inferredReadHp[readName] = 2;
                    inferredHp2Reads++;
                    readVoteRows.push_back({
                        chrom, readName, hp1, hp2,
                        minorityFraction, "HP2"
                    });
                }
            }
            germlineVoteSeconds += elapsedSeconds(voteT0);
        }

        const auto voteOutputT0 = std::chrono::steady_clock::now();
        std::sort(
            readVoteRows.begin(),
            readVoteRows.end(),
            [](const ReadHaplotypeVoteRow& a,
               const ReadHaplotypeVoteRow& b) {
                if (a.chrom != b.chrom) return a.chrom < b.chrom;
                return a.read_name < b.read_name;
            }
        );
        const fs::path readVotePath =
            fs::path(outDir) / "read_haplotype_votes.tsv";
        std::ofstream readVoteOut(readVotePath);
        if (!readVoteOut.is_open()) {
            std::cerr << "Error: cannot write " << readVotePath << "\n";
            return 1;
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
        germlineVoteSeconds += elapsedSeconds(voteOutputT0);
        std::cout << "[germline_haplotag] read_vote_tsv="
                  << readVotePath
                  << " rows=" << readVoteRows.size() << "\n";

        const std::size_t unassignedReads =
            candidateReadCount -
            inferredHp1Reads -
            inferredHp2Reads -
            mixedReads;
        const bool hasUsableAssignments =
            inferredHp1Reads + inferredHp2Reads > 0;
        std::cout << "[germline_haplotag]"
                  << " candidate_reads=" << candidateReadCount
                  << " min_informative_germline_mutations="
                  << kMinInformativeGermlineMutations
                  << " informative_observations="
                  << informativeGermlineObservations
                  << " hp1_reads=" << inferredHp1Reads
                  << " hp2_reads=" << inferredHp2Reads
                  << " mixed_reads=" << mixedReads
                  << " unassigned_reads=" << unassignedReads
                  << " has_hp_assignments="
                  << (hasUsableAssignments ? 1 : 0) << "\n";
        std::cout << "[haplotag_status] mode=germline_vcf"
                  << " haplotagged="
                  << (hasUsableAssignments ? 1 : 0)
                  << " source=germline-vcf\n";
    }

    // Second pass: attach BAM or inferred read HP values and run the existing
    // somatic annotation and pair-comparison logic.
    for (auto& chunk : chunks) {
        const auto parseT0 = std::chrono::steady_clock::now();
        const std::unordered_map<std::string, int>* externalReadHp = nullptr;
        if (useGermlineVcf) {
            auto chromHpIt = inferredReadHpByChrom.find(chunk.chrom);
            if (chromHpIt != inferredReadHpByChrom.end()) {
                externalReadHp = &chromHpIt->second;
            }
        }
        auto pileup = readMpileup(
            chunk.mpileup_out.string(),
            chunk.halo_snvs,
            !useGermlineVcf && bamIsHaplotagged,
            externalReadHp
        );
        somaticParseAttachSeconds += elapsedSeconds(parseT0);

        for (auto& snv : chunk.halo_snvs) {
            if (!snv.needs_read_vaf) continue;

            const std::string key =
                snv.chrom + ":" + std::to_string(snv.pos);
            auto cachedVaf = readVafByKey.find(key);
            if (cachedVaf != readVafByKey.end()) {
                snv.vaf = cachedVaf->second;
                continue;
            }
            if (skippedReadVafKeys.count(key) != 0) {
                pileup.erase(key);
                continue;
            }

            auto pileupIt = pileup.find(key);
            float readVaf = 0.0f;
            if (pileupIt != pileup.end() &&
                calculateVafFromReads(snv, pileupIt->second, readVaf)) {
                snv.vaf = readVaf;
                readVafByKey.emplace(key, readVaf);
            } else {
                skippedReadVafKeys.insert(key);
                hapSummaries.erase(key);
                pileup.erase(key);
                std::cerr << "Warning: skipping somatic mutation "
                          << key
                          << " because no pileup reads support REF or ALT\n";
            }
        }

        const auto summaryT0 = std::chrono::steady_clock::now();
        for (auto& snv : chunk.halo_snvs) {
            const std::string key =
                snv.chrom + ":" + std::to_string(snv.pos);
            auto pileupIt = pileup.find(key);
            if (pileupIt == pileup.end()) continue;

            annotateSNVHpSummary(snv, pileupIt->second);
            if (snv.hp1_alt_reads + snv.hp2_alt_reads <
                minHaplotaggedAltReadsPerSomaticMutation) {
                snv.hp_label = "UNKNOWN";
                snv.hp_proportion = 0.0f;
            }
            auto& summary = hapSummaries[key];
            summary.total_coverage = snv.total_coverage_reads;
            summary.hp1_ref_reads = snv.hp1_ref_reads;
            summary.hp1_alt_reads = snv.hp1_alt_reads;
            summary.hp2_ref_reads = snv.hp2_ref_reads;
            summary.hp2_alt_reads = snv.hp2_alt_reads;
            summary.nohp_alt_reads = snv.nohp_alt_reads;
        }
        somaticHpSummarySeconds += elapsedSeconds(summaryT0);

        const auto pairT0 = std::chrono::steady_clock::now();
        {
            std::ofstream chunkOut(chunk.chunk_out);
            std::size_t ownedPairs = 0;
            std::size_t comparedPairs = 0;
            std::size_t crossBoundaryPairs = 0;

            for (std::size_t i = 0; i < chunk.halo_snvs.size(); ++i) {
                if (chunk.halo_snvs[i].pos < chunk.region_start ||
                    chunk.halo_snvs[i].pos > chunk.region_end) {
                    continue;
                }

                for (std::size_t j = i + 1;
                     j < chunk.halo_snvs.size();
                     ++j) {
                    if (chunk.halo_snvs[j].pos -
                            chunk.halo_snvs[i].pos >
                        maxPairDistanceBp) {
                        break;
                    }
                    ownedPairs++;
                    if (chunk.halo_snvs[j].pos > chunk.region_end) {
                        crossBoundaryPairs++;
                    }

                    const std::string key1 =
                        chunk.halo_snvs[i].chrom + ":" +
                        std::to_string(chunk.halo_snvs[i].pos);
                    const std::string key2 =
                        chunk.halo_snvs[j].chrom + ":" +
                        std::to_string(chunk.halo_snvs[j].pos);
                    if (pileup.count(key1) && pileup.count(key2)) {
                        compareSNVs(
                            chunk.halo_snvs[i],
                            chunk.halo_snvs[j],
                            pileup[key1],
                            pileup[key2],
                            chunkOut
                        );
                        comparedPairs++;
                    }
                }
            }

            std::cout << "[chunk_halo_pairs] chrom=" << chunk.chrom
                      << " core_start=" << chunk.region_start
                      << " core_end=" << chunk.region_end
                      << " owned_distance_eligible=" << ownedPairs
                      << " compared=" << comparedPairs
                      << " cross_boundary_distance_eligible="
                      << crossBoundaryPairs << "\n";
        }
        processFile(
            chunk.chunk_out.string(),
            minReads,
            divergentSameHp
        );
        pairComparisonSeconds += elapsedSeconds(pairT0);
    }

    std::cout << "\n All processing complete.\n";
    std::cout << "   Pileup files in: " << pileupDir << "\n";
    std::cout << "   Chunk files in:  " << chunkDir << "\n";
    writeHaplotaggedSnvsTsv(
        fs::path(outDir),
        snvs,
        hapSummaries,
        skippedReadVafKeys,
        minHaplotaggedAltReadsPerSomaticMutation
    );

    int retainedSnvCount = 0;
    for (const auto& snv : snvs) {
        const std::string key =
            snv.chrom + ":" + std::to_string(snv.pos);
        if (skippedReadVafKeys.count(key) == 0) {
            retainedSnvCount++;
        }
    }
    std::cout << "[read_vaf_fallback] candidates="
              << readVafCandidateKeys.size()
              << " calculated=" << readVafByKey.size()
              << " skipped=" << skippedReadVafKeys.size() << "\n";

    const auto graphT0 = std::chrono::steady_clock::now();
    generateGraphs(outDir, retainedSnvCount);
    const double graphSeconds = elapsedSeconds(graphT0);

    const std::string timingMode =
        useGermlineVcf ? "germline_vcf" : "bam_hp_tags";
    std::cout << "[haplotag_timing] mode=" << timingMode
              << " stage=somatic_pileup seconds="
              << somaticPileupSeconds << "\n";
    std::cout << "[haplotag_timing] mode=" << timingMode
              << " stage=candidate_collection seconds="
              << candidateCollectionSeconds << "\n";
    std::cout << "[haplotag_timing] mode=" << timingMode
              << " stage=germline_pileup seconds="
              << germlinePileupSeconds << "\n";
    std::cout << "[haplotag_timing] mode=" << timingMode
              << " stage=germline_vote_resolution seconds="
              << germlineVoteSeconds << "\n";
    std::cout << "[haplotag_timing] mode=" << timingMode
              << " stage=somatic_parse_attach seconds="
              << somaticParseAttachSeconds << "\n";
    std::cout << "[haplotag_timing] mode=" << timingMode
              << " stage=somatic_hp_summary seconds="
              << somaticHpSummarySeconds << "\n";
    std::cout << "[haplotag_timing] mode=" << timingMode
              << " stage=pair_comparison seconds="
              << pairComparisonSeconds << "\n";
    std::cout << "[haplotag_timing] mode=" << timingMode
              << " stage=graph_generation seconds="
              << graphSeconds << "\n";
    std::cout << "[haplotag_timing] mode=" << timingMode
              << " stage=main_total seconds="
              << elapsedSeconds(main_t0) << "\n";

    return 0;
}
