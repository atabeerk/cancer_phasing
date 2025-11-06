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
        bcf_unpack(rec, BCF_UN_STR);

        SNV s;
        s.chrom = bcf_hdr_id2name(hdr, rec->rid);
        s.pos = rec->pos + 1;  // convert from 0-based to 1-based

        // Extract REF and ALT alleles (assumes single ALT)
        const char* ref = rec->d.allele[0];
        const char* alt = rec->d.allele[1];
        s.ref = ref[0];
        s.alt = alt[0];

        snvs.push_back(s);
    }

    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(fp);

    return snvs;
}


// ---- Parse mpileup (-s), but use VCF ref bases ----
std::unordered_map<std::string, PositionInfo> readMpileup(
    const std::string& mpileupFile,
    const std::vector<SNV>& snvs)
{
    // Create lookup: chrom:pos → (ref, alt)
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
        std::string chrom, bases, quality, readNamesStr;
        int pos, depth;
        char ref; // will be ignored (always 'N')

        iss >> chrom >> pos >> ref >> depth >> bases >> quality;
        std::getline(iss, readNamesStr);
        if (readNamesStr.empty()) continue;

        while (!readNamesStr.empty() && std::isspace(readNamesStr[0]))
            readNamesStr.erase(0, 1);

        // Split read names by comma
        std::vector<std::string> readNames;
        size_t start = 0, end = 0;
        while ((end = readNamesStr.find(',', start)) != std::string::npos) {
            readNames.push_back(readNamesStr.substr(start, end - start));
            start = end + 1;
        }
        readNames.push_back(readNamesStr.substr(start));

        // Parse the bases column — skip over +/− sequences
        std::string parsedBases;
        for (size_t i = 0; i < bases.size();) {
            char c = bases[i];

            if (c == '+' || c == '-') {
                // Parse the number following +/-
                size_t j = i + 1;
                std::string num;
                while (j < bases.size() && std::isdigit(bases[j])) {
                    num += bases[j++];
                }
                int len = std::stoi(num);
                // Skip the inserted/deleted bases
                j += len;
                i = j;
                continue;
            }

            // Ignore read start (^) and end ($) markers
            if (c == '^') { i += 2; continue; }
            if (c == '$') { i += 1; continue; }

            // Otherwise, it's a real base call
            parsedBases += c;
            ++i;
        }

        if (readNames.size() != parsedBases.size()) {
            std::cerr << "Warning: read count != parsed base count at "
                      << chrom << ":" << pos
                      << " (" << readNames.size() << " vs "
                      << parsedBases.size() << ")\n";
        }

        std::string key = chrom + ":" + std::to_string(pos);
        PositionInfo info;
        info.ref = refAltLookup[key].first;
        info.alt = refAltLookup[key].second;

        for (size_t i = 0; i < std::min(readNames.size(), parsedBases.size()); ++i)
            info.readBase[readNames[i]] = parsedBases[i];

        mp[key] = info;
    }

    return mp;
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

    out << s1.chrom << "\t" << s1.pos << "\t" << s2.pos
        << "\tALT_ALT=" << altAlt
        << "\tALT_REF=" << altRef
        << "\tREF_ALT=" << refAlt
        << "\tREF_REF=" << refRef
        << "\tTOTAL=" << (altAlt + altRef + refAlt + refRef) << "\n";
}
