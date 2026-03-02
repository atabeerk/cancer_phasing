#include "snv_analysis.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cctype>
#include <htslib/vcf.h>
#include <htslib/hts.h>

// ---- Parse VCF (.vcf or .vcf.gz) using htslib ----
std::vector<SNV> readVCF(const std::string& vcfFile) {
    std::vector<SNV> snvs;

    htsFile* fp = bcf_open(vcfFile.c_str(), "r");
    if (!fp) {
        std::cerr << "Error opening VCF file: " << vcfFile << "\n";
        return snvs;
    }

    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    if (!hdr) {
        std::cerr << "Error reading VCF header\n";
        bcf_close(fp);
        return snvs;
    }

    bcf1_t* rec = bcf_init();
    if (!rec) {
        std::cerr << "Error allocating VCF record\n";
        bcf_hdr_destroy(hdr);
        bcf_close(fp);
        return snvs;
    }

    while (bcf_read(fp, hdr, rec) == 0) {
        bcf_unpack(rec, BCF_UN_ALL);

        SNV s;
        s.chrom = bcf_hdr_id2name(hdr, rec->rid);
        s.pos = rec->pos + 1;  // convert from 0-based to 1-based

        // Extract REF and ALT alleles (assumes single ALT)
        const char* ref = rec->d.allele[0];
        const char* alt = rec->d.allele[1];
        s.ref = ref[0];
        s.alt = alt[0];

        float* af = NULL;
        int naf = 0;
        if (bcf_get_info_float(hdr, rec, "AF", &af, &naf) > 0 && naf > 0) {
            s.vaf = af[0];
            free(af);
        } else {
            // Try FORMAT/AD (allele depth) to calculate VAF
            int32_t* ad = NULL;
            int nad = 0;
            if (bcf_get_format_int32(hdr, rec, "AD", &ad, &nad) > 0 && nad >= 2) {
                int ref_depth = ad[0];
                int alt_depth = ad[1];
                int total_depth = ref_depth + alt_depth;
                s.vaf = (total_depth > 0) ? (float)alt_depth / total_depth : 0.5f;
                free(ad);
            } else {
                // Fallback: use 0.5 if VAF not available
                s.vaf = 0.5f;
                std::cerr << "Warning: No VAF/AF/AD info for " << s.chrom << ":" << s.pos 
                          << ", using default 0.5\n";
            }
        }

        snvs.push_back(s);
    }

    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(fp);

    return snvs;
}


// ---- Parse mpileup and read HP tags directly from output-extra columns ----
std::unordered_map<std::string, PositionInfo> readMpileup(
    const std::string& mpileupFile,
    const std::vector<SNV>& snvs)
{
    // Create lookup: chrom:pos -> (ref, alt)
    std::unordered_map<std::string, std::pair<char, char>> refAltLookup;
    for (const auto& s : snvs) {
        std::string key = s.chrom + ":" + std::to_string(s.pos);
        refAltLookup[key] = {s.ref, s.alt};
    }

    std::unordered_map<std::string, PositionInfo> mp;
    std::ifstream in(mpileupFile);
    if (!in.is_open()) {
        std::cerr << "Error opening mpileup: " << mpileupFile << "\n";
        return mp;
    }

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;

        std::istringstream iss(line);
        std::string chrom, bases, quality, readNamesStr, hpTagsStr;
        int pos, depth;
        char ref;

        iss >> chrom >> pos >> ref >> depth >> bases >> quality;
        iss >> readNamesStr;
        iss >> hpTagsStr;
        if (readNamesStr.empty()) continue;

        std::vector<std::string> readNames;
        size_t start = 0, end = 0;
        while ((end = readNamesStr.find(',', start)) != std::string::npos) {
            readNames.push_back(readNamesStr.substr(start, end - start));
            start = end + 1;
        }
        readNames.push_back(readNamesStr.substr(start));

        std::vector<std::string> hpTags;
        if (!hpTagsStr.empty()) {
            start = 0;
            end = 0;
            while ((end = hpTagsStr.find(',', start)) != std::string::npos) {
                hpTags.push_back(hpTagsStr.substr(start, end - start));
                start = end + 1;
            }
            hpTags.push_back(hpTagsStr.substr(start));
        }

        // Parse the bases column and skip +/- indel blocks.
        std::string parsedBases;
        for (size_t i = 0; i < bases.size();) {
            char c = bases[i];

            if (c == '+' || c == '-') {
                size_t j = i + 1;
                std::string num;
                while (j < bases.size() && std::isdigit(bases[j])) {
                    num += bases[j++];
                }
                int len = std::stoi(num);
                j += len;
                i = j;
                continue;
            }

            if (c == '^') { i += 2; continue; }
            if (c == '$') { i += 1; continue; }

            parsedBases += c;
            ++i;
        }

        if (readNames.size() != parsedBases.size()) {
            std::cerr << "Warning: read count != parsed base count at "
                      << chrom << ":" << pos
                      << " (" << readNames.size() << " vs " << parsedBases.size() << ")\n";
        }

        std::string key = chrom + ":" + std::to_string(pos);
        PositionInfo info;
        info.ref = refAltLookup[key].first;
        info.alt = refAltLookup[key].second;

        for (size_t i = 0; i < std::min(readNames.size(), parsedBases.size()); ++i) {
            info.readBase[readNames[i]] = parsedBases[i];
            if (i < hpTags.size()) {
                const std::string& hp = hpTags[i];
                if (hp == "1") info.readHP[readNames[i]] = 1;
                else if (hp == "2") info.readHP[readNames[i]] = 2;
            }
        }

        mp[key] = info;
    }

    return mp;
}


static std::string majorityHaplotypeForAltSupportingReads(const PositionInfo& p) {
    int hp1 = 0;
    int hp2 = 0;
    for (const auto& [read, base] : p.readBase) {
        bool isAlt = (std::toupper(base) == std::toupper(p.alt));
        if (!isAlt) continue;
        auto itHP = p.readHP.find(read);
        if (itHP == p.readHP.end()) continue;
        if (itHP->second == 1) hp1++;
        else if (itHP->second == 2) hp2++;
    }

    if (hp1 == 0 && hp2 == 0) return "UNKNOWN";
    if (hp1 == hp2) return "MIXED";
    return (hp1 > hp2) ? "HP1" : "HP2";
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

    const std::string hap1 = majorityHaplotypeForAltSupportingReads(p1);
    const std::string hap2 = majorityHaplotypeForAltSupportingReads(p2);

    out << s1.chrom << "\t" << s1.pos << "\t" << s2.pos
        << "\tVAF1=" << s1.vaf << "\tVAF2=" << s2.vaf
        << "\tALT_ALT=" << altAlt
        << "\tALT_REF=" << altRef
        << "\tREF_ALT=" << refAlt
        << "\tREF_REF=" << refRef
        << "\tTOTAL=" << (altAlt + altRef + refAlt + refRef)
        << "\tHAP1=" << hap1
        << "\tHAP2=" << hap2
        << "\n";
}
