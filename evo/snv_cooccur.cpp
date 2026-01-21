#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <unordered_map>
#include <vector>
#include <filesystem>

#include "htslib/vcf.h"
#include "htslib/sam.h"

#include "snv_analysis.hpp"
#include "process_snv_output.hpp"

namespace fs = std::filesystem;

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


void generateGraphs(const std::string& outDir) {
    std::string chunkDir = fs::path(outDir) / "chunk_files";
    std::string graphDir = fs::path(outDir) / "graphs";
    std::string cmd = "python graph_ops/main.py " + chunkDir + " --outdir " + graphDir;

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
    // ← MODIFIED: Removed MIN_READ_SUPPORT parameter ←
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <VCF_FILE> <BAM_FILE> <OUTPUT_DIR>\n";
        return 1;
    }

    std::string vcfFile = argv[1];
    std::string bamFile = argv[2];
    std::string outDir = argv[3];

    // Main output directory
    if (!fs::exists(outDir)) {
        fs::create_directories(outDir);
    }

    // Subdirectories for organization
    fs::path pileupDir = fs::path(outDir) / "pileup_files";
    fs::path chunkDir  = fs::path(outDir) / "chunk_files";

    fs::create_directories(pileupDir);
    fs::create_directories(chunkDir);

    // Read SNVs from VCF
    auto snvs = readVCF(vcfFile);

    // Group SNVs by chromosome
    std::unordered_map<std::string, std::vector<SNV>> chrom_snvs;
    for (const auto& s : snvs)
        chrom_snvs[s.chrom].push_back(s);

    const int chunk_size = 10'000'000;

    for (const auto& [chrom, snv_list] : chrom_snvs) {
        if (snv_list.empty()) continue;

        int min_pos = snv_list.front().pos;
        int max_pos = snv_list.back().pos;

        for (int region_start = min_pos; region_start <= max_pos; region_start += chunk_size) {
            int region_end = region_start + chunk_size - 1;

            // Extract SNVs in this region
            std::vector<SNV> chunk_snvs;
            for (const auto& s : snv_list) {
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
            std::string cmd = "samtools mpileup -q 0 -Q 0 --output-QNAME -l " + pos_file.string() +
                              " " + bamFile + " -o " + mpileup_out.string();

            std::cout << "Running command:\n" << cmd << "\n";
            if (system(cmd.c_str()) != 0) {
                std::cerr << "Error: mpileup failed for region " << chrom
                          << ":" << region_start << "-" << region_end << "\n";
                continue;
            }

            // --- Read mpileup output ---
            auto pileup = readMpileup(mpileup_out, chunk_snvs);

            // --- Write chunk comparison file under chunk_files ---
            fs::path chunk_out_file = chunkDir / ("chunk_" + chrom + "_" +
                                                  std::to_string(region_start) + "_" +
                                                  std::to_string(region_end) + ".txt");

            {
                std::ofstream chunk_out(chunk_out_file);
                for (size_t i = 0; i < chunk_snvs.size(); ++i) {
                    for (size_t j = i + 1; j < chunk_snvs.size(); ++j) {
                        std::string key1 = chunk_snvs[i].chrom + ":" + std::to_string(chunk_snvs[i].pos);
                        std::string key2 = chunk_snvs[j].chrom + ":" + std::to_string(chunk_snvs[j].pos);

                        if (pileup.count(key1) && pileup.count(key2)) {
                            compareSNVs(chunk_snvs[i], chunk_snvs[j], pileup[key1], pileup[key2], chunk_out);
                        }
                    }
                }
            }
            
            // ← MODIFIED: Removed minReadSupport parameter ←
            processFile(chunk_out_file.string());
        }
    }

    std::cout << "\n All processing complete.\n";
    std::cout << "   Pileup files in: " << pileupDir << "\n";
    std::cout << "   Chunk files in:  " << chunkDir  << "\n";

    generateGraphs(outDir);

    return 0;
}