#include "snv_analysis.hpp"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <cctype>
#include <htslib/vcf.h>
#include <htslib/hts.h>

// ---- Parse VCF (.vcf or .vcf.gz) using htslib ----
bool readVCF(const std::string& vcfFile, std::vector<SNV>& snvs) {
    snvs.clear();

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

    bcf1_t* rec = bcf_init();
    if (!rec) {
        std::cerr << "Error allocating VCF record\n";
        bcf_hdr_destroy(hdr);
        bcf_close(fp);
        return false;
    }

    int status = 0;
    while ((status = bcf_read(fp, hdr, rec)) == 0) {
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
        snvs.push_back(s);
    }

    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(fp);
    if (status < -1) std::cerr << "Error reading prepared somatic VCF: " << vcfFile << "\n";
    return status >= -1;
}

bool readGermlineVCF(const std::string& vcfFile, std::vector<GermlineSNV>& snvs) {
    snvs.clear();
    htsFile* fp = bcf_open(vcfFile.c_str(), "r");
    if (!fp) {
        std::cerr << "Error opening germline VCF file: " << vcfFile << "\n";
        return false;
    }
    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    if (!hdr || bcf_hdr_nsamples(hdr) != 1) {
        std::cerr << "Error: prepared germline VCF must contain one sample\n";
        if (hdr) bcf_hdr_destroy(hdr);
        bcf_close(fp);
        return false;
    }
    bcf1_t* rec = bcf_init();
    if (!rec) {
        bcf_hdr_destroy(hdr);
        bcf_close(fp);
        return false;
    }
    int32_t* genotypes = nullptr;
    int genotypeCapacity = 0;
    bool valid = true;
    int status = 0;
    while ((status = bcf_read(fp, hdr, rec)) == 0) {
        bcf_unpack(rec, BCF_UN_ALL);
        if (rec->n_allele != 2 ||
            bcf_get_genotypes(hdr, rec, &genotypes, &genotypeCapacity) < 2) {
            std::cerr << "Error: invalid marker in prepared germline VCF\n";
            valid = false;
            break;
        }
        snvs.push_back({
            bcf_hdr_id2name(hdr, rec->rid), static_cast<int>(rec->pos + 1),
            rec->d.allele[0][0], rec->d.allele[1][0],
            bcf_gt_allele(genotypes[0]) == 1 ? 1 : 2
        });
    }
    free(genotypes);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(fp);
    if (status < -1) std::cerr << "Error reading prepared germline VCF\n";
    return valid && status >= -1;
}

struct ParsedPileupRow {
    std::string chrom;
    int pos = 0;
    std::vector<std::string> readNames;
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

    constexpr std::size_t blockWidth = 4;
    if (fields.size() < 3 + blockWidth ||
        (fields.size() - 3) % blockWidth != 0) {
        return false;
    }

    row.chrom = fields[0];
    try {
        std::size_t parsed = 0;
        row.pos = std::stoi(fields[1], &parsed);
        if (parsed != fields[1].size() || row.pos < 1) return false;
    } catch (const std::exception&) {
        return false;
    }
    row.readNames.clear();
    row.bases.clear();

    for (std::size_t blockStart = 3;
         blockStart < fields.size();
         blockStart += blockWidth) {
        int depth = 0;
        try {
            std::size_t parsed = 0;
            depth = std::stoi(fields[blockStart], &parsed);
            if (parsed != fields[blockStart].size() || depth < 0) return false;
        } catch (const std::exception&) {
            return false;
        }
        if (depth == 0) {
            continue;
        }

        const std::string& rawBases = fields[blockStart + 1];
        const auto blockReadNames =
            splitCommaSeparated(fields[blockStart + 3]);
        if (blockReadNames.size() != static_cast<std::size_t>(depth)) return false;
        row.readNames.insert(
            row.readNames.end(),
            blockReadNames.begin(),
            blockReadNames.end()
        );

        const std::size_t baseStart = row.bases.size();
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
                if (number.empty()) return false;
                try {
                    const int length = std::stoi(number);
                    if (length < 1 || static_cast<std::size_t>(length) > rawBases.size() - j) return false;
                    i = j + static_cast<std::size_t>(length);
                } catch (const std::exception&) {
                    return false;
                }
                continue;
            }
            if (c == '^') {
                if (i + 1 >= rawBases.size()) return false;
                i += 2;
                continue;
            }
            if (c == '$') {
                ++i;
                continue;
            }
            row.bases += c;
            ++i;
        }
        if (row.bases.size() - baseStart != static_cast<std::size_t>(depth)) return false;
    }
    return true;
}

std::unordered_map<std::string, HaplotypeSnvSummary> summarizeSomaticMpileup(
    const std::string& mpileupFile,
    const std::vector<SNV>& snvs,
    const std::unordered_map<std::string, int>* externalReadHp
) {
    std::unordered_map<std::string, std::pair<char, char>> refAltLookup;
    for (const auto& snv : snvs) {
        refAltLookup[snv.chrom + ":" + std::to_string(snv.pos)] =
            {snv.ref, snv.alt};
    }

    std::unordered_map<std::string, HaplotypeSnvSummary> result;
    std::ifstream input(mpileupFile);
    if (!input.is_open()) {
        throw std::runtime_error("Cannot open mpileup: " + mpileupFile);
    }

    std::string line;
    while (std::getline(input, line)) {
        if (line.empty()) continue;

        ParsedPileupRow row;
        if (!parsePileupLine(line, row))
            throw std::runtime_error("Malformed somatic mpileup row in " + mpileupFile);
        const std::string key =
            row.chrom + ":" + std::to_string(row.pos);
        auto alleleIt = refAltLookup.find(key);
        if (alleleIt == refAltLookup.end()) continue;

        PositionInfo info;
        info.ref = alleleIt->second.first;
        info.alt = alleleIt->second.second;
        const std::size_t observationCount =
            row.readNames.size();
        for (std::size_t i = 0; i < observationCount; ++i) {
            const std::string& readName = row.readNames[i];
            info.readBase[readName] = row.bases[i];
        }
        result[key] = summarizeSnvReads(info, externalReadHp);
    }
    if (input.bad()) throw std::runtime_error("Cannot finish reading mpileup: " + mpileupFile);
    return result;
}

void collectReadNamesFromMpileup(
    const std::string& mpileupFile,
    std::unordered_set<std::string>& readNames
) {
    std::ifstream input(mpileupFile);
    if (!input.is_open()) {
        throw std::runtime_error("Cannot open mpileup: " + mpileupFile);
    }

    std::string line;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        ParsedPileupRow row;
        if (!parsePileupLine(line, row))
            throw std::runtime_error("Malformed mpileup row in " + mpileupFile);
        const std::size_t observationCount =
            row.readNames.size();
        for (std::size_t i = 0; i < observationCount; ++i) {
            readNames.insert(row.readNames[i]);
        }
    }
    if (input.bad()) throw std::runtime_error("Cannot finish reading mpileup: " + mpileupFile);
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
        throw std::runtime_error("Cannot open germline mpileup: " + mpileupFile);
    }

    std::size_t informativeObservations = 0;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        ParsedPileupRow row;
        if (!parsePileupLine(line, row))
            throw std::runtime_error("Malformed mpileup row in " + mpileupFile);

        const std::string key =
            row.chrom + ":" + std::to_string(row.pos);
        auto markerIt = markerLookup.find(key);
        if (markerIt == markerLookup.end()) continue;
        const GermlineSNV& marker = *markerIt->second;

        const std::size_t observationCount =
            row.readNames.size();
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
    if (input.bad()) throw std::runtime_error("Cannot finish reading germline mpileup: " + mpileupFile);
    return informativeObservations;
}


HaplotypeSnvSummary summarizeSnvReads(
    const PositionInfo& p, const std::unordered_map<std::string, int>* readHaplotypes
) {
    HaplotypeSnvSummary summary;
    summary.total_coverage = static_cast<int>(p.readBase.size());
    for (const auto& [read, base] : p.readBase) {
        const bool isAlt = (std::toupper(static_cast<unsigned char>(base)) == std::toupper(static_cast<unsigned char>(p.alt)));
        const bool isRef = (base == '.' || base == ',' ||
                            std::toupper(static_cast<unsigned char>(base)) == std::toupper(static_cast<unsigned char>(p.ref)));
        if (!isAlt && !isRef) continue;
        int hp = 0;
        if (readHaplotypes != nullptr) {
            const auto assignment = readHaplotypes->find(read);
            if (assignment != readHaplotypes->end()) hp = assignment->second;
        }
        if (hp == 1) {
            if (isAlt) summary.hp1_alt_reads++;
            else summary.hp1_ref_reads++;
        } else if (hp == 2) {
            if (isAlt) summary.hp2_alt_reads++;
            else summary.hp2_ref_reads++;
        } else {
            if (isAlt) summary.nohp_alt_reads++;
            else summary.nohp_ref_reads++;
        }
    }
    return summary;
}

bool calculateHaplotypeVaf(
    int refReads,
    int altReads,
    float& vaf
) {
    const int alleleDepth = refReads + altReads;
    if (alleleDepth == 0) return false;

    vaf = static_cast<float>(altReads) /
          static_cast<float>(alleleDepth);
    return true;
}
