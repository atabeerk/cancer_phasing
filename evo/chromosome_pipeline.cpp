#include <algorithm>
#include <cstdlib>
#include <iostream>

#include "pipeline.hpp"

namespace fs = std::filesystem;

namespace {
static constexpr int kMinInformativeGermlineMutations = 2;

class ChromosomePipeline {
public:
    explicit ChromosomePipeline(const PipelineOptions& options) : options(options),
        useGermlineVcf(!options.germlineVcfFile.empty()),
        pileupDir(fs::path(options.outDir) / "pileup_files"),
        chunkDir(fs::path(options.outDir) / "chunk_files") {}

    bool loadInputs();
    bool createSomaticPileups();
    bool inferReadHaplotypes();
    void summarizeMutations();
    void filterAndWriteMutations();
    bool inferRelationships();
    bool createGraphs();

private:
    const PipelineOptions& options;
    const bool useGermlineVcf;
    const fs::path pileupDir;
    const fs::path chunkDir;
    std::string bamArguments;
    std::vector<SNV> snvs;
    std::vector<GermlineSNV> germlineSnvs;
    std::unordered_map<std::string, std::vector<SNV>> chromSnvs;
    std::unordered_map<std::string, HaplotypeSnvSummary> hapSummaries;
    std::unordered_map<std::string, std::unordered_set<std::string>> candidateReadsByChrom;
    std::unordered_map<std::string, std::unordered_map<std::string, int>> inferredReadHpByChrom;
    int retainedSnvCount = 0;
};

bool ChromosomePipeline::loadInputs() {
    for (const auto& bamFile : options.bamFiles) {
        bamArguments += " " + shellQuote(bamFile);
    }
    std::cout << "[haplotag_mode] mode="
              << (useGermlineVcf ? "germline_vcf" : "none")
              << "\n";
    std::cout << "[relationship_filter] divergent_same_hp="
              << (options.divergentSameHp ? 1 : 0) << "\n";

    if (!fs::exists(options.outDir)) {
        fs::create_directories(options.outDir);
    }
    fs::create_directories(pileupDir);
    fs::create_directories(chunkDir);

    if (!readVCF(options.vcfFile, snvs)) return false;
    if (!filterExcludedRegions(snvs, options.excludeRegionsBedFiles)) return false;
    if (useGermlineVcf && !readGermlineVCF(options.germlineVcfFile, germlineSnvs)) {
        return false;
    }

    hapSummaries.reserve(snvs.size());
    for (const auto& snv : snvs) {
        chromSnvs[snv.chrom].push_back(snv);
    }

    if (!useGermlineVcf) {
        std::cout << "[haplotag_status] mode=none haplotagged=0 source=none\n";
    }

    return true;
}

bool ChromosomePipeline::createSomaticPileups() {
    for (auto& [chrom, mutations] : chromSnvs) {
        std::sort(mutations.begin(), mutations.end(),
                  [](const SNV& a, const SNV& b) { return a.pos < b.pos; });
        const fs::path positionsFile = pileupDir / ("positions_" + chrom + ".txt");
        const fs::path mpileupOut = pileupDir / ("mpileup_" + chrom + ".out");
        writePileupPositions(positionsFile, mutations);
        const std::string region = chrom + ":" +
            std::to_string(mutations.front().pos) + "-" +
            std::to_string(mutations.back().pos);
        const std::string command = mpileupCommand(
            bamArguments, region, positionsFile, mpileupOut
        );

        std::cout << "[somatic_pileup] chrom=" << chrom
                  << " snvs=" << mutations.size() << "\n";
        std::cout << "Running command:\n" << command << "\n";

        const int pileupStatus = system(command.c_str());
        if (pileupStatus != 0) {
            std::cerr << "Error: mpileup failed for chromosome " << chrom << "\n";
            return false;
        }

        if (useGermlineVcf) {
            collectReadNamesFromMpileup(mpileupOut.string(), candidateReadsByChrom[chrom]);
        }
    }
    return true;
}

bool ChromosomePipeline::inferReadHaplotypes() {
    if (!useGermlineVcf) return true;
    std::size_t candidateReadCount = 0;
    std::size_t informativeGermlineObservations = 0;
    std::size_t inferredHp1Reads = 0;
    std::size_t inferredHp2Reads = 0;
    std::size_t mixedReads = 0;
    std::vector<ReadHaplotypeVoteRow> readVoteRows;

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
        writePileupPositions(positionsFile, markers);

        const fs::path mpileupOut =
            pileupDir / ("germline_mpileup_" + chrom + ".out");
        const std::string region = chrom + ":" +
            std::to_string(markers.front().pos) + "-" +
            std::to_string(markers.back().pos);
        const std::string command = mpileupCommand(
            bamArguments, region, positionsFile, mpileupOut
        );

        std::cout << "Running germline command:\n"
                  << command << "\n";

        const int pileupStatus = system(command.c_str());
        if (pileupStatus != 0) {
            std::cerr << "Error: germline mpileup failed for "
                      << chrom << "\n";
            return false;
        }

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
            const std::string assignment = hpAssignmentFromCounts(
                hp1, hp2, kMinInformativeGermlineMutations
            );
            if (assignment == "UNKNOWN") {
                readVoteRows.push_back({
                    chrom, readName, hp1, hp2, -1.0, "UNASSIGNED"
                });
                continue;
            }

            const int minor = std::min(hp1, hp2);
            const double minorityFraction =
                static_cast<double>(minor) /
                static_cast<double>(total);
            if (assignment == "MIXED") {
                mixedReads++;
            } else {
                const int hp = assignment == "HP1" ? 1 : 2;
                inferredReadHp[readName] = hp;
                if (hp == 1) inferredHp1Reads++;
                else inferredHp2Reads++;
            }
            readVoteRows.push_back({
                chrom, readName, hp1, hp2, minorityFraction, assignment
            });
        }
    }

    const fs::path readVotePath = fs::path(options.outDir) / "read_haplotype_votes.tsv";
    if (!writeReadHaplotypeVotes(readVotePath, readVoteRows)) return false;

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
    return true;
}

void ChromosomePipeline::summarizeMutations() {
    for (const auto& [chrom, mutations] : chromSnvs) {
        const std::unordered_map<std::string, int>* externalReadHp = nullptr;
        if (useGermlineVcf) {
            auto hp = inferredReadHpByChrom.find(chrom);
            if (hp != inferredReadHpByChrom.end()) externalReadHp = &hp->second;
        }
        const fs::path pileup = pileupDir / ("mpileup_" + chrom + ".out");
        auto summaries = summarizeSomaticMpileup(
            pileup.string(), mutations, externalReadHp
        );
        hapSummaries.merge(summaries);
    }
}

void ChromosomePipeline::filterAndWriteMutations() {
    // Filter using unique reads per mutation.
    // Use the same REF + ALT denominator as the read-derived mutation VAF.
    const std::size_t inputSupportSnvs = snvs.size();
    snvs.erase(
        std::remove_if(
            snvs.begin(), snvs.end(),
            [this](const SNV& snv) {
                const std::string key =
                    snv.chrom + ":" + std::to_string(snv.pos);
                const auto it = hapSummaries.find(key);
                if (it == hapSummaries.end()) return true;
                const auto& summary = it->second;
                const int altReads = summary.altReadCount();
                const int alleleDepth = summary.alleleDepth();
                // Strict OR: exactly 3 ALT reads or exactly 8% alone is
                // insufficient. Integer arithmetic preserves the 8% boundary.
                const bool keep = altReads > 3 ||
                    (alleleDepth > 0 && 100LL * altReads > 8LL * alleleDepth);
                if (!keep) hapSummaries.erase(it);
                return !keep;
            }
        ),
        snvs.end()
    );
    std::cout << "[somatic_support_filter] min_alt_reads_exclusive=3"
              << " min_alt_percent_exclusive=8 logic=or coverage=ref_plus_alt"
              << " input_snvs=" << inputSupportSnvs
              << " excluded_snvs=" << inputSupportSnvs - snvs.size()
              << " retained_snvs=" << snvs.size() << "\n";

    std::cout << "\n=== Writing mutation summaries ===\n";
    std::cout << "   Pileup files in: " << pileupDir << "\n";
    std::cout << "   Chunk files in:  " << chunkDir << "\n";
    writeSnvCountsTsv(
        fs::path(options.outDir),
        snvs,
        hapSummaries
    );

    retainedSnvCount = static_cast<int>(snvs.size());
}

bool ChromosomePipeline::inferRelationships() {
    const fs::path inferenceScript = fs::absolute(options.executablePath).parent_path() /
        "inference/infer_evo_relationships.py";
    std::string inferenceCommand = "python " + shellQuote(inferenceScript.string()) +
        " --vcf " + shellQuote(options.vcfFile) +
        " --output-dir " + shellQuote(options.outDir) +
        " --evo-pileup-bam-count " + std::to_string(options.bamFiles.size()) +
        " --window-size " + std::to_string(options.windowSize) +
        " --window-stride " + std::to_string(options.windowStride);
    if (options.divergentSameHp) {
        inferenceCommand += " --divergent-same-hp";
    }
    std::cout << "[relationship_inference] method=shared_tree_"
              << options.windowSize << "ms" << options.windowStride
              << " maximum_mutations=" << options.windowSize
              << " group_stride=" << options.windowStride
              << " minimum_informative_reads=3 minimum_alt_alt_reads=2"
              << " graph_support_gate=1 em_iterations=50\n";
    std::cout << "Running command: " << inferenceCommand << std::endl;

    if (std::system(inferenceCommand.c_str()) != 0) {
        std::cerr << "Error: shared-tree relationship inference failed\n";
        return false;
    }
    return true;
}

bool ChromosomePipeline::createGraphs() {
    return generateGraphs(options, retainedSnvCount);
}
} // namespace

int runChromosomePipeline(const PipelineOptions& options) {
    ChromosomePipeline pipeline(options);
    if (!pipeline.loadInputs()) return 1;
    if (!pipeline.createSomaticPileups()) return 1;
    if (!pipeline.inferReadHaplotypes()) return 1;
    pipeline.summarizeMutations();
    pipeline.filterAndWriteMutations();
    if (!pipeline.inferRelationships()) return 1;
    if (!pipeline.createGraphs()) return 1;
    if (!options.germlineVcfFile.empty()) {
        const fs::path outDir(options.outDir);
        fs::rename(outDir / "snv_counts.tsv", outDir / "haplotagged_snvs.tsv");
    }
    return 0;
}
