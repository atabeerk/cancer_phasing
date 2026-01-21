#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <map>
#include <vector>

using namespace std;

struct SNPEntry {
    string chr;
    long pos1;
    long pos2;
    int ALT_ALT = 0;
    int ALT_REF = 0;
    int REF_ALT = 0;
    int REF_REF = 0;
    int TOTAL = 0;
    bool valid = true;
    
    float vaf1 = 0.5f;
    float vaf2 = 0.5f;
};

SNPEntry parseLine(const string& line) {
    SNPEntry e;
    stringstream ss(line);
    if (!(ss >> e.chr >> e.pos1 >> e.pos2)) {
        e.valid = false;
        return e;
    }

    string kv;
    while (ss >> kv) {
        auto pos = kv.find('=');
        if (pos == string::npos) continue;
        string key = kv.substr(0, pos);
        string val_str = kv.substr(pos + 1);
        
        if (key == "VAF1") {
            e.vaf1 = stof(val_str);
        } else if (key == "VAF2") {
            e.vaf2 = stof(val_str);
        } else if (key == "ALT_ALT") {
            e.ALT_ALT = stoi(val_str);
        } else if (key == "ALT_REF") {
            e.ALT_REF = stoi(val_str);
        } else if (key == "REF_ALT") {
            e.REF_ALT = stoi(val_str);
        } else if (key == "REF_REF") {
            e.REF_REF = stoi(val_str);
        } else if (key == "TOTAL") {
            e.TOTAL = stoi(val_str);
        }
    }

    if (e.TOTAL <= 0) e.valid = false;
    return e;
}


float calculateReliabilityForRelationship(const SNPEntry& e,
                                          const std::string& relType,
                                          int minReads)
{
    if (e.TOTAL == 0) return 0.0f;

    const float obs_RR = static_cast<float>(e.REF_REF) / e.TOTAL;
    const float obs_AR = static_cast<float>(e.ALT_REF) / e.TOTAL;
    const float obs_RA = static_cast<float>(e.REF_ALT) / e.TOTAL;
    const float obs_AA = static_cast<float>(e.ALT_ALT) / e.TOTAL;

    float reliability = 0.0f;

    if (relType == "snp1_before_snp2") {
        // requires VAF1 >= VAF2
        if (e.vaf1 < e.vaf2) {
            return 0.0f;
        }

        // Required: R/R, A/R, A/A
        if (e.REF_REF < minReads || e.ALT_REF < minReads || e.ALT_ALT < minReads) {
            return 0.0f;
        }

        // Relationship-aware expected fractions:
        float exp_RR = 1.0f - e.vaf1;
        float exp_AR = e.vaf1 - e.vaf2;
        float exp_AA = e.vaf2;

        if (exp_RR <= 0.0f || exp_AR <= 0.0f || exp_AA <= 0.0f) {
            return 0.0f;
        }

        float min_ratio = 1.0f;
        min_ratio = std::min(min_ratio, obs_RR / exp_RR);
        min_ratio = std::min(min_ratio, obs_AR / exp_AR);
        min_ratio = std::min(min_ratio, obs_AA / exp_AA);

        // Forbidden: R/A
        float forbidden_penalty = obs_RA;

        reliability = min_ratio - forbidden_penalty;
        return std::max(0.0f, reliability);
    }

    if (relType == "snp2_before_snp1") {
        // requires VAF2 >= VAF1
        if (e.vaf2 < e.vaf1) {
            return 0.0f;
        }

        // Required: R/R, R/A, A/A
        if (e.REF_REF < minReads || e.REF_ALT < minReads || e.ALT_ALT < minReads) {
            return 0.0f;
        }

        // Relationship-aware expected fractions:
        float exp_RR = 1.0f - e.vaf2;
        float exp_RA = e.vaf2 - e.vaf1;
        float exp_AA = e.vaf1;

        if (exp_RR <= 0.0f || exp_RA <= 0.0f || exp_AA <= 0.0f) {
            return 0.0f;
        }

        float min_ratio = 1.0f;
        min_ratio = std::min(min_ratio, obs_RR / exp_RR);
        min_ratio = std::min(min_ratio, obs_RA / exp_RA);
        min_ratio = std::min(min_ratio, obs_AA / exp_AA);

        // Forbidden: A/R
        float forbidden_penalty = obs_AR;

        reliability = min_ratio - forbidden_penalty;
        return std::max(0.0f, reliability);
    }

    if (relType == "cooccurring") {
        // Required: R/R, A/A
        if (e.REF_REF < minReads || e.ALT_ALT < minReads) {
            return 0.0f;
        }

        float exp_RR = std::min(1.0f - e.vaf1, 1.0f - e.vaf2);
        float exp_AA = std::min(e.vaf1, e.vaf2);

        float min_ratio = 1.0f;
        if (exp_RR > 0.0f) min_ratio = std::min(min_ratio, obs_RR / exp_RR);
        if (exp_AA > 0.0f) min_ratio = std::min(min_ratio, obs_AA / exp_AA);

        // Forbidden: A/R, R/A
        float forbidden_penalty = obs_AR + obs_RA;

        reliability = min_ratio - forbidden_penalty;
        return std::max(0.0f, reliability);
    }

    // Co-occurring with loss
    if (relType == "cooccurring_loss") {
        // Required: A/A only
        if (e.ALT_ALT < minReads) {
            return 0.0f;
        }

        float exp_AA = std::min(e.vaf1, e.vaf2);

        float min_ratio = 1.0f;
        if (exp_AA > 0.0f) min_ratio = std::min(min_ratio, obs_AA / exp_AA);

        // Forbidden: R/R, A/R, R/A
        float forbidden_penalty = obs_RR + obs_AR + obs_RA;

        reliability = min_ratio - forbidden_penalty;
        return std::max(0.0f, reliability);
    }

    // Divergent
    if (relType == "divergent") {
        // Required: A/R, R/A
        if (e.ALT_REF < minReads || e.REF_ALT < minReads) {
            return 0.0f;
        }

        float exp_AR = std::min(e.vaf1, 1.0f - e.vaf2);
        float exp_RA = std::min(1.0f - e.vaf1, e.vaf2);

        float min_ratio = 1.0f;
        if (exp_AR > 0.0f) min_ratio = std::min(min_ratio, obs_AR / exp_AR);
        if (exp_RA > 0.0f) min_ratio = std::min(min_ratio, obs_RA / exp_RA);

        // Forbidden: R/R, A/A
        float forbidden_penalty = obs_RR + obs_AA;

        reliability = min_ratio - forbidden_penalty;
        return std::max(0.0f, reliability);
    }

    // SNP1 before SNP2 with loss (SNP1=A branch only)
    if (relType == "snp1_before_snp2_loss") {
        // Feasibility: SNP1-before-SNP2-loss requires VAF1 >= VAF2
        if (e.vaf1 < e.vaf2) {
            return 0.0f;
        }

        // Required: A/R, A/A
        if (e.ALT_REF < minReads || e.ALT_ALT < minReads) {
            return 0.0f;
        }

        // Expected fractions in SNP1=A branch:
        float exp_AR = e.vaf1 - e.vaf2;
        float exp_AA = e.vaf2;

        if (exp_AR <= 0.0f || exp_AA <= 0.0f) {
            return 0.0f;
        }

        float min_ratio = 1.0f;
        min_ratio = std::min(min_ratio, obs_AR / exp_AR);
        min_ratio = std::min(min_ratio, obs_AA / exp_AA);

        // Forbidden: R/R, R/A
        float forbidden_penalty = obs_RR + obs_RA;

        reliability = min_ratio - forbidden_penalty;
        return std::max(0.0f, reliability);
    }

    // SNP2 before SNP1 with loss (SNP2=A branch only)
    if (relType == "snp2_before_snp1_loss") {
        // requires VAF2 >= VAF1
        if (e.vaf2 < e.vaf1) {
            return 0.0f;
        }

        // Required: R/A, A/A
        if (e.REF_ALT < minReads || e.ALT_ALT < minReads) {
            return 0.0f;
        }

        float exp_RA = e.vaf2 - e.vaf1;
        float exp_AA = e.vaf1;

        if (exp_RA <= 0.0f || exp_AA <= 0.0f) {
            return 0.0f;
        }

        float min_ratio = 1.0f;
        min_ratio = std::min(min_ratio, obs_RA / exp_RA);
        min_ratio = std::min(min_ratio, obs_AA / exp_AA);

        // Forbidden: R/R, A/R
        float forbidden_penalty = obs_RR + obs_AR;

        reliability = min_ratio - forbidden_penalty;
        return std::max(0.0f, reliability);
    }

    // Unknown relation type
    return 0.0f;
}

struct Classification {
    string type;
    float final_reliability;
    float best_score;
    float margin;
};


// Number of reads that "support" a given relation type
int supportingReadsForRelation(const SNPEntry& e, const std::string& relType) {
    if (relType == "snp1_before_snp2") {
        // Required: R/R, A/R, A/A
        return e.REF_REF + e.ALT_REF + e.ALT_ALT;
    }
    if (relType == "snp2_before_snp1") {
        // Required: R/R, R/A, A/A
        return e.REF_REF + e.REF_ALT + e.ALT_ALT;
    }
    if (relType == "cooccurring") {
        // Required: R/R, A/A
        return e.REF_REF + e.ALT_ALT;
    }
    if (relType == "cooccurring_loss") {
        // Required: A/A only
        return e.ALT_ALT;
    }
    if (relType == "divergent") {
        // Required: A/R, R/A
        return e.ALT_REF + e.REF_ALT;
    }
    if (relType == "snp1_before_snp2_loss") {
        // Required: A/R, A/A
        return e.ALT_REF + e.ALT_ALT;
    }
    if (relType == "snp2_before_snp1_loss") {
        // Required: R/A, A/A
        return e.REF_ALT + e.ALT_ALT;
    }
    return 0;
}


// Linear coverage weight: w = s / (s + SCALE)
float coverageWeight(int support, float scale) {
    if (support <= 0 || scale <= 0.0f) {
        return 0.0f;
    }
    float s = static_cast<float>(support);
    return s / (s + scale);
}


Classification classifyEntry(const SNPEntry& e, int minReads, float coverageScale) {
    // Calculate reliability for all relationship types (pattern-based)
    map<string, float> reliabilities;
    reliabilities["snp1_before_snp2"] = calculateReliabilityForRelationship(e, "snp1_before_snp2", minReads);
    reliabilities["snp2_before_snp1"] = calculateReliabilityForRelationship(e, "snp2_before_snp1", minReads);
    reliabilities["cooccurring"] = calculateReliabilityForRelationship(e, "cooccurring", minReads);
    reliabilities["cooccurring_loss"] = calculateReliabilityForRelationship(e, "cooccurring_loss", minReads);
    reliabilities["divergent"] = calculateReliabilityForRelationship(e, "divergent", minReads);
    reliabilities["snp1_before_snp2_loss"] = calculateReliabilityForRelationship(e, "snp1_before_snp2_loss", minReads);
    reliabilities["snp2_before_snp1_loss"]  = calculateReliabilityForRelationship(e, "snp2_before_snp1_loss", minReads);

    // Find best and second-best pattern scores
    string best_type = "error";
    float best_score = -999.0f;
    float second_best = -999.0f;

    for (const auto& kv : reliabilities) {
        const string& type  = kv.first;
        float score         = kv.second;
        if (score > best_score) {
            second_best = best_score;
            best_score  = score;
            best_type   = type;
        } else if (score > second_best) {
            second_best = score;
        }
    }

    // Calculate margin and base reliability from pattern matching
    float margin = best_score - second_best;
    float base_reliability = best_score * (1.0f + margin);
    float final_reliability = base_reliability;

    // If best score is <= 0, mark as error
    if (best_score <= 0.0f) {
        best_type = "error";
        final_reliability = 0.0f;
    } else {
        // Weight by the number of reads that support the chosen relation,
        // scaled by the median coverage (coverageScale) for this file.
        int support  = supportingReadsForRelation(e, best_type);
        float weight = coverageWeight(support, coverageScale);
        final_reliability = base_reliability * weight;
    }

    Classification result;
    result.type = best_type;
    result.final_reliability = final_reliability;
    result.best_score = best_score;
    result.margin = margin;

    return result;
}


void processFile(const string& filename, int minReads) {
    // First pass: read all entries to estimate a typical coverage (median TOTAL)
    ifstream infile(filename);
    if (!infile) {
        cerr << "Error: cannot open input file.\n";
        return;
    }

    vector<int> coverages;
    string line;

    while (getline(infile, line)) {
        if (line.empty()) continue;

        SNPEntry e = parseLine(line);
        if (!e.valid) continue;

        // Use TOTAL read count as a coverage proxy
        coverages.push_back(e.TOTAL);
    }

    infile.close();

    if (coverages.empty()) {
        cerr << "Warning: no valid entries found in " << filename << "\n";
        return;
    }

    // Compute median coverage to use as SCALE in weight = s / (s + SCALE)
    sort(coverages.begin(), coverages.end());
    int medianIndex = static_cast<int>(coverages.size() / 2);
    float coverageScale = static_cast<float>(coverages[medianIndex]);

    // Second pass: classify entries and write outputs using coverageScale
    infile.open(filename);
    if (!infile) {
        cerr << "Error: cannot reopen input file.\n";
        return;
    }

    string base = filename.substr(0, filename.find_last_of('.'));
    ofstream f_co(base + "_cooccurring.txt");
    ofstream f_div(base + "_divergent.txt");
    ofstream f_snp1(base + "_snp1_before_snp2.txt");
    ofstream f_snp2(base + "_snp2_before_snp1.txt");
    ofstream f_snp1_loss(base + "_snp1_before_snp2_loss.txt");
    ofstream f_snp2_loss(base + "_snp2_before_snp1_loss.txt");
    ofstream f_co_loss(base + "_cooccurring_loss.txt");
    ofstream f_err(base + "_errors.txt");

    while (getline(infile, line)) {
        if (line.empty()) continue;

        SNPEntry e = parseLine(line);
        if (!e.valid) continue;

        Classification cls = classifyEntry(e, minReads, coverageScale);

        // Format output with reliability scores
        string formatted =
            e.chr + string(" ") + to_string(e.pos1) + " " + to_string(e.pos2) + " "
            "VAF1=" + to_string(e.vaf1) + " "
            "VAF2=" + to_string(e.vaf2) + " "
            "ALT_ALT=" + to_string(e.ALT_ALT) + " "
            "ALT_REF=" + to_string(e.ALT_REF) + " "
            "REF_ALT=" + to_string(e.REF_ALT) + " "
            "REF_REF=" + to_string(e.REF_REF) + " "
            "TOTAL=" + to_string(e.TOTAL) + " "
            "RELIABILITY=" + to_string(cls.final_reliability) + " "
            "BEST_SCORE=" + to_string(cls.best_score) + " "
            "MARGIN=" + to_string(cls.margin);

        // Route to appropriate file based on classification
        if (cls.type == "cooccurring") {
            f_co << formatted << '\n';
        } else if (cls.type == "divergent") {
            f_div << formatted << '\n';
        } else if (cls.type == "snp1_before_snp2") {
            f_snp1 << formatted << '\n';
        } else if (cls.type == "snp2_before_snp1") {
            f_snp2 << formatted << '\n';
        } else if (cls.type == "snp1_before_snp2_loss") {
            f_snp1_loss << formatted << '\n';
        } else if (cls.type == "snp2_before_snp1_loss") {
            f_snp2_loss << formatted << '\n';
        } else if (cls.type == "cooccurring_loss") {
            f_co_loss << formatted << '\n';
        } else {
            // cls.type == "error" or final_reliability <= 0
            f_err << formatted << '\n';
        }
    }

    infile.close();
    f_co.close();
    f_div.close();
    f_snp1.close();
    f_snp2.close();
    f_snp1_loss.close();
    f_snp2_loss.close();
    f_co_loss.close();
    f_err.close();

    cout << "Processing complete.\n";
}