#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <unordered_map>
#include <vector>
#include <filesystem>
#include <algorithm>
#include <chrono>
#include <random>
#include <cstdint>

#include "htslib/vcf.h"
#include "htslib/sam.h"

#include "snv_analysis.hpp"
#include "process_snv_output.hpp"

namespace fs = std::filesystem;


bool commandAvailable(const std::string& cmd) {
    const std::string probe = "command -v " + cmd + " >/dev/null 2>&1";
    return std::system(probe.c_str()) == 0;
}

int detectMaxPairDistanceFromBam(const std::string& bamFile) {
    const auto t0 = std::chrono::steady_clock::now();
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

    int max_ref_span = 0;
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

bool detectBamHaplotagged(const std::string& bamFile, int sampleReads = 1000) {
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

    std::vector<uint8_t> sampled_has_hp;
    sampled_has_hp.reserve(sampleReads > 0 ? sampleReads : 0);
    std::uint64_t total_seen = 0;
    std::mt19937_64 rng(std::random_device{}());

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

    int hp_tagged_sampled_reads = 0;
    for (uint8_t v : sampled_has_hp) {
        if (v != 0) hp_tagged_sampled_reads++;
    }
    const bool has_hp_tag = (hp_tagged_sampled_reads > 0);

    bam_destroy1(rec);
    bam_hdr_destroy(hdr);
    sam_close(fp);

    std::cout << "[haplotag_detect] sampled_reads=" << sampled_has_hp.size()
              << " total_seen=" << total_seen
              << " hp_tagged_sampled_reads=" << hp_tagged_sampled_reads
              << " has_hp_tag=" << (has_hp_tag ? 1 : 0) << "\n";
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


int main(int argc, char* argv[]) {
    auto printUsage = [&]() {
        std::cerr << "Usage: " << argv[0]
                  << " --vcf <VCF_FILE> --bam <BAM_FILE> --output-dir <OUTPUT_DIR>"
                  << " [--min-reads <N>] [--vcf-sample-name <SAMPLE_NAME>]"
                  << " [--max-pair-distance-kb <N>] [--bam-haplotagged <0|1>]\n";
    };

    std::string vcfFile;
    std::string bamFile;
    std::string outDir;
    int minReads = 2;
    std::string vcfSampleName;
    int maxPairDistanceKbArg = -1; // <=0 means auto-detect from BAM
    int bamHaplotaggedArg = -1;    // -1 auto-detect, 0 no, 1 yes

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            printUsage();
            return 0;
        } else if (arg == "--vcf" && i + 1 < argc) {
            vcfFile = argv[++i];
        } else if (arg == "--bam" && i + 1 < argc) {
            bamFile = argv[++i];
        } else if (arg == "--output-dir" && i + 1 < argc) {
            outDir = argv[++i];
        } else if (arg == "--min-reads" && i + 1 < argc) {
            minReads = std::stoi(argv[++i]);
        } else if (arg == "--vcf-sample-name" && i + 1 < argc) {
            vcfSampleName = argv[++i];
        } else if (arg == "--max-pair-distance-kb" && i + 1 < argc) {
            maxPairDistanceKbArg = std::stoi(argv[++i]);
            if (maxPairDistanceKbArg <= 0) {
                std::cerr << "Error: --max-pair-distance-kb must be > 0\n";
                return 1;
            }
        } else if (arg == "--bam-haplotagged" && i + 1 < argc) {
            bamHaplotaggedArg = std::stoi(argv[++i]);
            if (bamHaplotaggedArg != 0 && bamHaplotaggedArg != 1) {
                std::cerr << "Error: --bam-haplotagged must be 0 or 1\n";
                return 1;
            }
        } else if (arg.rfind("--", 0) == 0) {
            std::cerr << "Error: unknown or incomplete argument: " << arg << "\n";
            printUsage();
            return 1;
        } else {
            std::cerr << "Error: positional argument not supported: " << arg << "\n";
            printUsage();
            return 1;
        }
    }

    if (vcfFile.empty() || bamFile.empty() || outDir.empty()) {
        std::cerr << "Error: missing required arguments.\n";
        printUsage();
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
    if (!fs::exists(bamFile)) {
        std::cerr << "Error: BAM file not found: " << bamFile << "\n";
        return 1;
    }

    std::cout << "Using MIN_READS = " << minReads << "\n";

    // Main output directory
    if (!fs::exists(outDir)) {
        fs::create_directories(outDir);
    }

    // Subdirectories for organization
    fs::path pileupDir = fs::path(outDir) / "pileup_files";
    fs::path chunkDir  = fs::path(outDir) / "chunk_files";

    fs::create_directories(pileupDir);
    fs::create_directories(chunkDir);

    // Read SNVs from VCF (optionally using FORMAT/AD from a specific sample)
    auto snvs = readVCF(vcfFile, vcfSampleName);

    // Group SNVs by chromosome
    std::unordered_map<std::string, std::vector<SNV>> chrom_snvs;
    for (const auto& s : snvs)
        chrom_snvs[s.chrom].push_back(s);

    int max_pair_distance_bp = 0;
    if (maxPairDistanceKbArg > 0) {
        max_pair_distance_bp = maxPairDistanceKbArg * 1000;
        std::cout << "Using max SNV pair distance = " << max_pair_distance_bp
                  << " bp (from --max-pair-distance-kb)\n";
    } else {
        max_pair_distance_bp = detectMaxPairDistanceFromBam(bamFile);
        std::cout << "Using max SNV pair distance = " << max_pair_distance_bp
                  << " bp (rounded from max aligned reference span in BAM)\n";
    }

    bool bam_is_haplotagged = false;
    if (bamHaplotaggedArg == 0 || bamHaplotaggedArg == 1) {
        bam_is_haplotagged = (bamHaplotaggedArg == 1);
        std::cout << "Using BAM haplotagged status = " << (bam_is_haplotagged ? 1 : 0)
                  << " (from --bam-haplotagged)\n";
    } else {
        bam_is_haplotagged = detectBamHaplotagged(bamFile);
        std::cout << "Using BAM haplotagged status = " << (bam_is_haplotagged ? 1 : 0)
                  << " (auto-detected)\n";
    }

    const int chunk_size = 10'000'000;

    for (const auto& [chrom, snv_list] : chrom_snvs) {
        if (snv_list.empty()) continue;
        std::vector<SNV> sorted_snv_list = snv_list;
        std::sort(sorted_snv_list.begin(), sorted_snv_list.end(),
                  [](const SNV& a, const SNV& b) { return a.pos < b.pos; });

        int min_pos = sorted_snv_list.front().pos;
        int max_pos = sorted_snv_list.back().pos;

        for (int region_start = min_pos; region_start <= max_pos; region_start += chunk_size) {
            int region_end = region_start + chunk_size - 1;

            // Extract SNVs in this region
            std::vector<SNV> chunk_snvs;
            for (const auto& s : sorted_snv_list) {
                if (s.pos >= region_start && s.pos <= region_end)
                    chunk_snvs.push_back(s);
            }
            if (chunk_snvs.empty()) continue;

            // --- Create BED/position file under pileup_files ---
            fs::path pos_file = pileupDir / ("positions_" + chrom + "_" + std::to_string(region_start) + ".txt");
            {
                std::ofstream bed(pos_file);
                for (const auto& s : chunk_snvs)
                    bed << s.chrom << "\t" << (s.pos - 1) << "\t" << s.pos << "\n";
            }

            // --- Run mpileup and store its output under pileup_files ---
            fs::path mpileup_out = pileupDir / ("mpileup_" + chrom + "_" + std::to_string(region_start) + ".out");
            std::string cmd = "samtools mpileup -q 0 -Q 0 --output-QNAME";
            if (bam_is_haplotagged) {
                cmd += " --output-extra HP";
            }
            cmd += " -l " + pos_file.string() +
                              " " + bamFile + " -o " + mpileup_out.string();

            std::cout << "Running command:\n" << cmd << "\n";
            if (system(cmd.c_str()) != 0) {
                std::cerr << "Error: mpileup failed for region " << chrom
                          << ":" << region_start << "-" << region_end << "\n";
                continue;
            }

            // --- Read mpileup output ---
            auto pileup = readMpileup(mpileup_out, chunk_snvs);

            // --- Precompute per-SNV HP labels/read counts once per chunk ---
            for (auto& s : chunk_snvs) {
                std::string key = s.chrom + ":" + std::to_string(s.pos);
                auto it = pileup.find(key);
                if (it != pileup.end()) {
                    annotateSNVHpSummary(s, it->second);
                }
            }

            // --- Write chunk comparison file under chunk_files ---
            fs::path chunk_out_file = chunkDir / ("chunk_" + chrom + "_" +
                                                  std::to_string(region_start) + "_" +
                                                  std::to_string(region_end) + ".txt");

            {
                std::ofstream chunk_out(chunk_out_file);
                for (size_t i = 0; i < chunk_snvs.size(); ++i) {
                    for (size_t j = i + 1; j < chunk_snvs.size(); ++j) {
                        if (chunk_snvs[j].pos - chunk_snvs[i].pos > max_pair_distance_bp) {
                            break;
                        }
                        std::string key1 = chunk_snvs[i].chrom + ":" + std::to_string(chunk_snvs[i].pos);
                        std::string key2 = chunk_snvs[j].chrom + ":" + std::to_string(chunk_snvs[j].pos);

                        if (pileup.count(key1) && pileup.count(key2)) {
                            compareSNVs(chunk_snvs[i], chunk_snvs[j], pileup[key1], pileup[key2], chunk_out);
                        }
                    }
                }
            }
            
            processFile(chunk_out_file.string(), minReads);
        }
    }

    std::cout << "\n All processing complete.\n";
    std::cout << "   Pileup files in: " << pileupDir << "\n";
    std::cout << "   Chunk files in:  " << chunkDir  << "\n";

    generateGraphs(outDir, static_cast<int>(snvs.size()));

    return 0;
}