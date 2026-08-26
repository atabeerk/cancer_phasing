#include "snv_analysis.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cctype>
#include <algorithm>
#include <cstring>
#include <htslib/vcf.h>
#include <htslib/hts.h>

// ---- Parse VCF (.vcf or .vcf.gz) using htslib ----
static bool resolveSample(
    bcf_hdr_t* hdr,
    const std::string& requestedSampleName,
    bool requireExplicit,
    const std::string& label,
    int& selectedSampleIdx,
    std::string& resolvedSampleName
) {
    const int numSamples = bcf_hdr_nsamples(hdr);
    selectedSampleIdx = -1;
    resolvedSampleName.clear();

    if (requireExplicit && requestedSampleName.empty()) {
        std::cerr << "Error: " << label << " sample name is required.\n";
        return false;
    }

    if (!requestedSampleName.empty()) {
        selectedSampleIdx =
            bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, requestedSampleName.c_str());
        if (selectedSampleIdx < 0) {
            std::cerr << "Error: " << label << " sample '"
                      << requestedSampleName << "' was not found in the VCF.\n";
            return false;
        }
    } else if (numSamples == 1) {
        selectedSampleIdx = 0;
    } else if (numSamples > 1) {
        std::cerr << "Error: " << label << " VCF contains " << numSamples
                  << " samples; an explicit sample name is required.\n";
        return false;
    }

    if (selectedSampleIdx >= 0) {
        resolvedSampleName = hdr->samples[selectedSampleIdx];
        std::cout << "Using " << label << " VCF sample '"
                  << resolvedSampleName << "'.\n";
    }
    return true;
}

bool readVCF(
    const std::string& vcfFile,
    const std::string& requestedSampleName,
    std::vector<SNV>& snvs,
    std::string& resolvedSampleName
) {
    snvs.clear();
    resolvedSampleName.clear();

    htsFile* fp = bcf_open(vcfFile.c_str(), "r");
    if (!fp) {
        std::cerr << "Error opening VCF file: " << vcfFile << "\n";
        return false;
    }

    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    if (!hdr) {
        std::cerr << "Error reading VCF header\n";
        bcf_close(fp);
        return false;
    }

    int selectedSampleIdx = -1;
    if (!resolveSample(
            hdr,
            requestedSampleName,
            false,
            "somatic",
            selectedSampleIdx,
            resolvedSampleName)) {
        bcf_hdr_destroy(hdr);
        bcf_close(fp);
        return false;
    }

    bcf1_t* rec = bcf_init();
    if (!rec) {
        std::cerr << "Error allocating VCF record\n";
        bcf_hdr_destroy(hdr);
        bcf_close(fp);
        return false;
    }

    const int numSamples = bcf_hdr_nsamples(hdr);
    while (bcf_read(fp, hdr, rec) == 0) {
        bcf_unpack(rec, BCF_UN_ALL);
        if (rec->n_allele < 2 || rec->d.allele[0][0] == '\0' ||
            rec->d.allele[1][0] == '\0') {
            continue;
        }

        SNV s;
        s.chrom = bcf_hdr_id2name(hdr, rec->rid);
        s.pos = rec->pos + 1;
        s.ref = rec->d.allele[0][0];
        s.alt = rec->d.allele[1][0];

        bool vafSet = false;
        if (selectedSampleIdx >= 0) {
            int32_t* ad = nullptr;
            int nad = 0;
            const int nvalsTotal =
                bcf_get_format_int32(hdr, rec, "AD", &ad, &nad);
            if (nvalsTotal > 0 && numSamples > 0) {
                const int valsPerSample = nvalsTotal / numSamples;
                const int idx0 = selectedSampleIdx * valsPerSample;
                if (valsPerSample >= 2 && idx0 + 1 < nvalsTotal) {
                    const int32_t refDepth = ad[idx0];
                    const int32_t altDepth = ad[idx0 + 1];
                    const bool refOk =
                        refDepth != bcf_int32_missing &&
                        refDepth != bcf_int32_vector_end;
                    const bool altOk =
                        altDepth != bcf_int32_missing &&
                        altDepth != bcf_int32_vector_end;
                    if (refOk && altOk) {
                        const int totalDepth =
                            static_cast<int>(refDepth + altDepth);
                        if (totalDepth > 0) {
                            s.vaf = static_cast<float>(altDepth) /
                                    static_cast<float>(totalDepth);
                            vafSet = true;
                        }
                    }
                }
            }
            if (ad != nullptr) free(ad);
        }

        if (!vafSet) {
            float* af = nullptr;
            int naf = 0;
            if (bcf_get_info_float(hdr, rec, "AF", &af, &naf) > 0 &&
                naf > 0) {
                s.vaf = af[0];
                vafSet = true;
            }
            if (af != nullptr) free(af);
        }

        if (!vafSet) {
            s.needs_read_vaf = true;
            std::cerr << "Warning: No usable VAF/AF/AD at "
                      << s.chrom << ":" << s.pos
                      << ", will calculate VAF from somatic pileup reads\n";
        }
        snvs.push_back(s);
    }

    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(fp);
    return true;
}

bool readGermlineVCF(
    const std::string& vcfFile,
    const std::string& sampleName,
    std::vector<GermlineSNV>& snvs,
    GermlineVcfStats& stats
) {
    snvs.clear();
    stats = GermlineVcfStats{};

    htsFile* fp = bcf_open(vcfFile.c_str(), "r");
    if (!fp) {
        std::cerr << "Error opening germline VCF file: " << vcfFile << "\n";
        return false;
    }

    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    if (!hdr) {
        std::cerr << "Error reading germline VCF header\n";
        bcf_close(fp);
        return false;
    }

    int sampleIdx = -1;
    std::string resolvedSampleName;
    if (!resolveSample(
            hdr,
            sampleName,
            true,
            "germline",
            sampleIdx,
            resolvedSampleName)) {
        bcf_hdr_destroy(hdr);
        bcf_close(fp);
        return false;
    }

    bcf1_t* rec = bcf_init();
    if (!rec) {
        std::cerr << "Error allocating germline VCF record\n";
        bcf_hdr_destroy(hdr);
        bcf_close(fp);
        return false;
    }

    struct CandidateMarker {
        GermlineSNV snv;
        bool isPass;
    };
    std::vector<CandidateMarker> candidates;
    int32_t* genotypes = nullptr;
    int genotypeCapacity = 0;
    const int numSamples = bcf_hdr_nsamples(hdr);

    while (bcf_read(fp, hdr, rec) == 0) {
        bcf_unpack(rec, BCF_UN_ALL);
        stats.total_records++;

        const bool isPass =
            bcf_has_filter(hdr, rec, const_cast<char*>("PASS")) == 1;
        if (isPass) stats.pass_records++;

        if (rec->n_allele != 2 ||
            std::strlen(rec->d.allele[0]) != 1 ||
            std::strlen(rec->d.allele[1]) != 1) {
            continue;
        }

        const int genotypeValues = bcf_get_genotypes(
            hdr, rec, &genotypes, &genotypeCapacity);
        if (genotypeValues <= 0 || numSamples <= 0) continue;
        const int ploidy = genotypeValues / numSamples;
        if (ploidy < 2) continue;

        int32_t* sampleGt = genotypes + sampleIdx * ploidy;
        if (bcf_gt_is_missing(sampleGt[0]) ||
            bcf_gt_is_missing(sampleGt[1]) ||
            !bcf_gt_is_phased(sampleGt[1])) {
            continue;
        }

        const int allele0 = bcf_gt_allele(sampleGt[0]);
        const int allele1 = bcf_gt_allele(sampleGt[1]);
        int altHaplotype = 0;
        if (allele0 == 1 && allele1 == 0) {
            altHaplotype = 1;
        } else if (allele0 == 0 && allele1 == 1) {
            altHaplotype = 2;
        } else {
            continue;
        }

        candidates.push_back({
            GermlineSNV{
                bcf_hdr_id2name(hdr, rec->rid),
                static_cast<int>(rec->pos + 1),
                rec->d.allele[0][0],
                rec->d.allele[1][0],
                altHaplotype
            },
            isPass
        });
    }

    if (genotypes != nullptr) free(genotypes);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(fp);

    const bool requirePass = stats.pass_records > 0;
    snvs.reserve(candidates.size());
    for (const auto& candidate : candidates) {
        if (!requirePass || candidate.isPass) {
            snvs.push_back(candidate.snv);
        }
    }
    stats.eligible_records = static_cast<int>(snvs.size());

    std::cout << "[germline_vcf] total_records=" << stats.total_records
              << " pass_records=" << stats.pass_records
              << " require_pass=" << (requirePass ? 1 : 0)
              << " eligible_records=" << stats.eligible_records
              << " skipped_records="
              << (stats.total_records - stats.eligible_records) << "\n";
    return true;
}

struct ParsedPileupRow {
    std::string chrom;
    int pos = 0;
    std::vector<std::string> readNames;
    std::vector<std::string> hpTags;
    std::string bases;
};

static std::vector<std::string> splitCommaSeparated(
    const std::string& value
) {
    std::vector<std::string> fields;
    if (value.empty()) return fields;

    std::size_t start = 0;
    while (true) {
        const std::size_t end = value.find(',', start);
        if (end == std::string::npos) {
            fields.push_back(value.substr(start));
            break;
        }
        fields.push_back(value.substr(start, end - start));
        start = end + 1;
    }
    return fields;
}

static bool parsePileupLine(
    const std::string& line,
    bool parseHpTags,
    ParsedPileupRow& row
) {
    std::vector<std::string> fields;
    std::size_t fieldStart = 0;
    while (true) {
        const std::size_t fieldEnd = line.find('\t', fieldStart);
        if (fieldEnd == std::string::npos) {
            fields.push_back(line.substr(fieldStart));
            break;
        }
        fields.push_back(line.substr(fieldStart, fieldEnd - fieldStart));
        fieldStart = fieldEnd + 1;
    }

    const std::size_t blockWidth = parseHpTags ? 5 : 4;
    if (fields.size() < 3 + blockWidth ||
        (fields.size() - 3) % blockWidth != 0) {
        return false;
    }

    row.chrom = fields[0];
    try {
        row.pos = std::stoi(fields[1]);
    } catch (const std::exception&) {
        return false;
    }
    row.readNames.clear();
    row.hpTags.clear();
    row.bases.clear();

    for (std::size_t blockStart = 3;
         blockStart < fields.size();
         blockStart += blockWidth) {
        int depth = 0;
        try {
            depth = std::stoi(fields[blockStart]);
        } catch (const std::exception&) {
            return false;
        }
        if (depth == 0) {
            continue;
        }

        const std::string& rawBases = fields[blockStart + 1];
        const auto blockReadNames =
            splitCommaSeparated(fields[blockStart + 3]);
        row.readNames.insert(
            row.readNames.end(),
            blockReadNames.begin(),
            blockReadNames.end()
        );
        if (parseHpTags) {
            const auto blockHpTags =
                splitCommaSeparated(fields[blockStart + 4]);
            row.hpTags.insert(
                row.hpTags.end(),
                blockHpTags.begin(),
                blockHpTags.end()
            );
        }

        for (std::size_t i = 0; i < rawBases.size();) {
            const char c = rawBases[i];
            if (c == '+' || c == '-') {
                std::size_t j = i + 1;
                std::string number;
                while (j < rawBases.size() &&
                       std::isdigit(
                           static_cast<unsigned char>(rawBases[j])
                       )) {
                    number += rawBases[j++];
                }
                if (number.empty()) {
                    ++i;
                    continue;
                }
                const int length = std::stoi(number);
                i = std::min(
                    rawBases.size(),
                    j + static_cast<std::size_t>(length)
                );
                continue;
            }
            if (c == '^') {
                i = std::min(rawBases.size(), i + 2);
                continue;
            }
            if (c == '$') {
                ++i;
                continue;
            }
            row.bases += c;
            ++i;
        }
    }
    return true;
}

std::unordered_map<std::string, PositionInfo> readMpileup(
    const std::string& mpileupFile,
    const std::vector<SNV>& snvs,
    bool parseBamHpTags,
    const std::unordered_map<std::string, int>* externalReadHp
) {
    std::unordered_map<std::string, std::pair<char, char>> refAltLookup;
    for (const auto& snv : snvs) {
        refAltLookup[snv.chrom + ":" + std::to_string(snv.pos)] =
            {snv.ref, snv.alt};
    }

    std::unordered_map<std::string, PositionInfo> result;
    std::ifstream input(mpileupFile);
    if (!input.is_open()) {
        std::cerr << "Error opening mpileup: " << mpileupFile << "\n";
        return result;
    }

    std::string line;
    while (std::getline(input, line)) {
        if (line.empty()) continue;

        ParsedPileupRow row;
        if (!parsePileupLine(line, parseBamHpTags, row)) continue;
        const std::string key =
            row.chrom + ":" + std::to_string(row.pos);
        auto alleleIt = refAltLookup.find(key);
        if (alleleIt == refAltLookup.end()) continue;

        if (row.readNames.size() != row.bases.size()) {
            std::cerr << "Warning: read count != parsed base count at "
                      << row.chrom << ":" << row.pos << " ("
                      << row.readNames.size() << " vs "
                      << row.bases.size() << ")\n";
        }

        PositionInfo info;
        info.ref = alleleIt->second.first;
        info.alt = alleleIt->second.second;
        const std::size_t observationCount =
            std::min(row.readNames.size(), row.bases.size());
        for (std::size_t i = 0; i < observationCount; ++i) {
            const std::string& readName = row.readNames[i];
            info.readBase[readName] = row.bases[i];

            if (externalReadHp != nullptr) {
                auto hpIt = externalReadHp->find(readName);
                if (hpIt != externalReadHp->end() &&
                    (hpIt->second == 1 || hpIt->second == 2)) {
                    info.readHP[readName] = hpIt->second;
                }
            } else if (parseBamHpTags && i < row.hpTags.size()) {
                if (row.hpTags[i] == "1") {
                    info.readHP[readName] = 1;
                } else if (row.hpTags[i] == "2") {
                    info.readHP[readName] = 2;
                }
            }
        }
        result[key] = std::move(info);
    }
    return result;
}

void collectReadNamesFromMpileup(
    const std::string& mpileupFile,
    std::unordered_set<std::string>& readNames
) {
    std::ifstream input(mpileupFile);
    if (!input.is_open()) {
        std::cerr << "Error opening mpileup: " << mpileupFile << "\n";
        return;
    }

    std::string line;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        ParsedPileupRow row;
        if (!parsePileupLine(line, false, row)) continue;
        const std::size_t observationCount =
            std::min(row.readNames.size(), row.bases.size());
        for (std::size_t i = 0; i < observationCount; ++i) {
            readNames.insert(row.readNames[i]);
        }
    }
}

std::size_t accumulateGermlineVotes(
    const std::string& mpileupFile,
    const std::vector<GermlineSNV>& snvs,
    const std::unordered_set<std::string>& candidateReads,
    std::unordered_map<std::string, GermlineVote>& votes
) {
    std::unordered_map<std::string, const GermlineSNV*> markerLookup;
    markerLookup.reserve(snvs.size());
    for (const auto& snv : snvs) {
        markerLookup[snv.chrom + ":" + std::to_string(snv.pos)] = &snv;
    }

    std::ifstream input(mpileupFile);
    if (!input.is_open()) {
        std::cerr << "Error opening germline mpileup: "
                  << mpileupFile << "\n";
        return 0;
    }

    std::size_t informativeObservations = 0;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        ParsedPileupRow row;
        if (!parsePileupLine(line, false, row)) continue;

        const std::string key =
            row.chrom + ":" + std::to_string(row.pos);
        auto markerIt = markerLookup.find(key);
        if (markerIt == markerLookup.end()) continue;
        const GermlineSNV& marker = *markerIt->second;

        const std::size_t observationCount =
            std::min(row.readNames.size(), row.bases.size());
        for (std::size_t i = 0; i < observationCount; ++i) {
            const std::string& readName = row.readNames[i];
            if (candidateReads.find(readName) == candidateReads.end()) {
                continue;
            }

            const char base = row.bases[i];
            const bool supportsAlt =
                std::toupper(static_cast<unsigned char>(base)) ==
                std::toupper(static_cast<unsigned char>(marker.alt));
            const bool supportsRef =
                base == '.' || base == ',' ||
                std::toupper(static_cast<unsigned char>(base)) ==
                std::toupper(static_cast<unsigned char>(marker.ref));
            if (!supportsAlt && !supportsRef) continue;

            const int haplotype = supportsAlt
                ? marker.alt_haplotype
                : (marker.alt_haplotype == 1 ? 2 : 1);
            GermlineVote& vote = votes[readName];
            if (haplotype == 1) {
                vote.hp1++;
            } else {
                vote.hp2++;
            }
            informativeObservations++;
        }
    }
    return informativeObservations;
}


static void hpSummaryForAltSupportingReads(
    const PositionInfo& p,
    std::string& hp_label,
    float& hp_proportion,
    int& hp1_ref_reads,
    int& hp1_alt_reads,
    int& hp2_ref_reads,
    int& hp2_alt_reads,
    int& nohp_alt_reads,
    int& total_coverage_reads
) {
    int hp1_ref = 0;
    int hp1 = 0;
    int hp2_ref = 0;
    int hp2 = 0;
    int nohp = 0;
    total_coverage_reads = static_cast<int>(p.readBase.size());
    for (const auto& [read, base] : p.readBase) {
        bool isAlt = (std::toupper(base) == std::toupper(p.alt));
        bool isRef = (base == '.' || base == ',' ||
                      std::toupper(base) == std::toupper(p.ref));
        if (!isAlt && !isRef) continue;
        auto itHP = p.readHP.find(read);
        if (itHP == p.readHP.end()) {
            if (isAlt) nohp++;
            continue;
        }
        if (itHP->second == 1) {
            if (isAlt) hp1++;
            else hp1_ref++;
        } else if (itHP->second == 2) {
            if (isAlt) hp2++;
            else hp2_ref++;
        } else if (isAlt) {
            nohp++;
        }
    }

    hp1_ref_reads = hp1_ref;
    hp1_alt_reads = hp1;
    hp2_ref_reads = hp2_ref;
    hp2_alt_reads = hp2;
    nohp_alt_reads = nohp;
    const int total = hp1 + hp2;
    if (total == 0) {
        hp_label = "UNKNOWN";
        hp_proportion = 0.0f;
        return;
    }
    const int minor = std::min(hp1, hp2);
    if (100LL * minor > 15LL * total) {
        hp_label = "MIXED";
        hp_proportion = static_cast<float>(std::max(hp1, hp2)) /
                        static_cast<float>(total);
        return;
    }
    hp_label = (hp1 > hp2) ? "HP1" : "HP2";
    hp_proportion = static_cast<float>(std::max(hp1, hp2)) / static_cast<float>(total);
}

void annotateSNVHpSummary(SNV& snv, const PositionInfo& p) {
    hpSummaryForAltSupportingReads(
        p,
        snv.hp_label,
        snv.hp_proportion,
        snv.hp1_ref_reads,
        snv.hp1_alt_reads,
        snv.hp2_ref_reads,
        snv.hp2_alt_reads,
        snv.nohp_alt_reads,
        snv.total_coverage_reads
    );
}

bool calculateVafFromReads(
    const SNV& snv,
    const PositionInfo& p,
    float& vaf
) {
    int refReads = 0;
    int altReads = 0;
    for (const auto& [read, base] : p.readBase) {
        (void)read;
        const bool isAlt =
            std::toupper(base) == std::toupper(snv.alt);
        const bool isRef =
            base == '.' || base == ',' ||
            std::toupper(base) == std::toupper(snv.ref);
        if (isAlt) {
            altReads++;
        } else if (isRef) {
            refReads++;
        }
    }

    const int alleleDepth = refReads + altReads;
    if (alleleDepth == 0) return false;

    vaf = static_cast<float>(altReads) /
          static_cast<float>(alleleDepth);
    return true;
}

void compareSNVs(const SNV& s1, const SNV& s2,
                 const PositionInfo& p1, const PositionInfo& p2,
                 std::ostream& out)
{
    int altAlt = 0;
    int altRef = 0;
    int refAlt = 0;
    int refRef = 0;

    for (const auto& [read, base1] : p1.readBase) {
        auto it2 = p2.readBase.find(read);

        if (it2 == p2.readBase.end()) {
            continue;
        }

        char base2 = it2->second;

        bool isAlt1 = (std::toupper(base1) == std::toupper(p1.alt));
        bool isRef1 = (std::toupper(base1) == std::toupper(p1.ref));
        bool isAlt2 = (std::toupper(base2) == std::toupper(p2.alt));
        bool isRef2 = (std::toupper(base2) == std::toupper(p2.ref));

        // skip ambiguous bases
        if ((!isAlt1 && !isRef1) || (!isAlt2 && !isRef2))
           continue;

        if (isAlt1 && isAlt2) altAlt++;
        else if (isAlt1 && !isAlt2) altRef++;
        else if (!isAlt1 && isAlt2) refAlt++;
        else if (!isAlt1 && !isAlt2) refRef++;
    }

    std::string hap1 = s1.hp_label;
    std::string hap2 = s2.hp_label;
    int hp1_reads1 = s1.hp1_alt_reads;
    int hp2_reads1 = s1.hp2_alt_reads;
    int hp1_reads2 = s2.hp1_alt_reads;
    int hp2_reads2 = s2.hp2_alt_reads;

    out << s1.chrom << "\t" << s1.pos << "\t" << s2.pos
        << "\tVAF1=" << s1.vaf << "\tVAF2=" << s2.vaf
        << "\tALT_ALT=" << altAlt
        << "\tALT_REF=" << altRef
        << "\tREF_ALT=" << refAlt
        << "\tREF_REF=" << refRef
        << "\tTOTAL=" << (altAlt + altRef + refAlt + refRef)
        << "\tHAP1=" << hap1
        << "\tHAP2=" << hap2
        << "\tHP_READS1=" << hp1_reads1 << "/" << hp2_reads1
        << "\tHP_READS2=" << hp1_reads2 << "/" << hp2_reads2
        << "\n";
}
