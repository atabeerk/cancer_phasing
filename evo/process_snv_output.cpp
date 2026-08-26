#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <map>
#include <vector>
#include <cmath>
#include <algorithm>

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
    string hap1 = "UNKNOWN";
    string hap2 = "UNKNOWN";
    string hp_reads1 = "0/0";
    string hp_reads2 = "0/0";
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
        } else if (key == "HAP1") {
            e.hap1 = val_str;
        } else if (key == "HAP2") {
            e.hap2 = val_str;
        } else if (key == "HP_READS1") {
            e.hp_reads1 = val_str;
        } else if (key == "HP_READS2") {
            e.hp_reads2 = val_str;
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
    const int informative_reads = e.ALT_ALT + e.ALT_REF + e.REF_ALT;
    if (informative_reads == 0) return 0.0f;

    // Condition every relation on reads carrying at least one ALT allele.
    // REF/REF reads are uninformative about the relationship between the two
    // mutations and are excluded from both support and penalties.
    const float obs_AR = static_cast<float>(e.ALT_REF) / informative_reads;
    const float obs_RA = static_cast<float>(e.REF_ALT) / informative_reads;
    const float obs_AA = static_cast<float>(e.ALT_ALT) / informative_reads;

    if (relType == "snp1_before_snp2") {
        if (e.vaf1 < e.vaf2) return 0.0f;
        if (e.ALT_REF < minReads || e.ALT_ALT < minReads) return 0.0f;

        float exp_AR = e.vaf1 - e.vaf2;
        float exp_AA = e.vaf2;
        const float expected_total = exp_AR + exp_AA;
        if (expected_total <= 0.0f) return 0.0f;
        exp_AR /= expected_total;
        exp_AA /= expected_total;

        float min_ratio = 1.0f;
        min_ratio = std::min(min_ratio, obs_AR / exp_AR);
        min_ratio = std::min(min_ratio, obs_AA / exp_AA);
        return std::max(0.0f, min_ratio - obs_RA);
    }

    if (relType == "snp2_before_snp1") {
        if (e.vaf2 < e.vaf1) return 0.0f;
        if (e.REF_ALT < minReads || e.ALT_ALT < minReads) return 0.0f;

        float exp_RA = e.vaf2 - e.vaf1;
        float exp_AA = e.vaf1;
        const float expected_total = exp_RA + exp_AA;
        if (expected_total <= 0.0f) return 0.0f;
        exp_RA /= expected_total;
        exp_AA /= expected_total;

        float min_ratio = 1.0f;
        min_ratio = std::min(min_ratio, obs_RA / exp_RA);
        min_ratio = std::min(min_ratio, obs_AA / exp_AA);
        return std::max(0.0f, min_ratio - obs_AR);
    }

    if (relType == "cooccurring") {
        if (e.ALT_ALT < minReads) return 0.0f;
        return std::max(0.0f, obs_AA - obs_AR - obs_RA);
    }

    if (relType == "divergent") {
        if (e.ALT_REF < minReads || e.REF_ALT < minReads) return 0.0f;

        float exp_AR = std::min(e.vaf1, 1.0f - e.vaf2);
        float exp_RA = std::min(1.0f - e.vaf1, e.vaf2);
        const float expected_total = exp_AR + exp_RA;
        if (expected_total <= 0.0f) return 0.0f;
        exp_AR /= expected_total;
        exp_RA /= expected_total;

        float min_ratio = 1.0f;
        min_ratio = std::min(min_ratio, obs_AR / exp_AR);
        min_ratio = std::min(min_ratio, obs_RA / exp_RA);
        return std::max(0.0f, min_ratio - obs_AA);
    }

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
    if (relType == "snp1_before_snp2") return e.ALT_REF + e.ALT_ALT;
    if (relType == "snp2_before_snp1") return e.REF_ALT + e.ALT_ALT;
    if (relType == "cooccurring") return e.ALT_ALT;
    if (relType == "divergent") return e.ALT_REF + e.REF_ALT;
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


static inline float clamp01(float x) {
    return (x < 0.0f) ? 0.0f : ((x > 1.0f) ? 1.0f : x);
}


// Symmetric ratio r >= 1, with a floor on the denominator to avoid blow-ups at tiny VAF.
static inline float vaf_ratio(float v1, float v2) {
    const float vmin_floor = 0.02f;   // denominator floor for ratio
    
    float lo = std::min(v1, v2);
    float hi = std::max(v1, v2);
    lo = std::max(lo, vmin_floor);
    return hi / lo; // always >= 1
}

// Piecewise linear: 1 until r0, then linearly to 0 at rmax.
static inline float cooccurring_prior_from_ratio(float r) {
    const float r0   = 1.0f;   // coocc penalty if vafs are different
    const float rmax = 1.35f;   // coocc prior -> 0 by 35% ratio diff

    if (r <= r0) return 1.0f;
    if (r >= rmax) return 0.0f;
    return 1.0f - (r - r0) / (rmax - r0);
}


Classification classifyEntry(const SNPEntry& e, int minReads, float coverageScale) {
    // Calculate reliability for all relationship types (pattern-based)
    map<string, float> reliabilities;
    reliabilities["snp1_before_snp2"] = calculateReliabilityForRelationship(e, "snp1_before_snp2", minReads);
    reliabilities["snp2_before_snp1"] = calculateReliabilityForRelationship(e, "snp2_before_snp1", minReads);
    reliabilities["cooccurring"] = calculateReliabilityForRelationship(e, "cooccurring", minReads);
    reliabilities["divergent"] = calculateReliabilityForRelationship(e, "divergent", minReads);

    // --- Compute ratio-based distance ---
    const float r = vaf_ratio(e.vaf1, e.vaf2);

    // --- Compute priors ---
    const float eps  = 0.05f;   // timing prior floor (prevents hard zero)

    const float p_co = clamp01(cooccurring_prior_from_ratio(r));
    const float p_tim = clamp01(eps + (1.0f - eps) * (1.0f - p_co));
    const float prior_strength = 1.0f; // 0=no prior effect, 1=full prior effect
    const float p_co_eff = 1.0f - prior_strength * (1.0f - p_co);
    const float p_tim_eff = 1.0f - prior_strength * (1.0f - p_tim);

    // --- Apply priors (multiplicative) ---
    // Cooccurring gets boosted when VAFs are similar (p_co close to 1)
    reliabilities["cooccurring"]      *= p_co_eff;

    // Timing gets punished when VAFs are similar
    reliabilities["snp1_before_snp2"]      *= p_tim_eff;
    reliabilities["snp2_before_snp1"]      *= p_tim_eff;

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
        // Weight by ALT-informative support, scaled by the median ALT-informative
        // coverage for this file.
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


static bool hasSameAssignedHaplotype(const SNPEntry& e) {
    return (e.hap1 == "HP1" && e.hap2 == "HP1") ||
           (e.hap1 == "HP2" && e.hap2 == "HP2");
}


void processFile(
    const string& filename,
    int minReads,
    bool divergentSameHp
) {
    // First pass: read all entries to estimate a typical coverage (median TOTAL)
    ifstream infile(filename);
    if (!infile) {
        cerr << "Error: cannot open input file.\n";
        return;
    }

    vector<int> coverages;
    string line;

    int retainedDivergent = 0;
    int filteredDivergent = 0;

    while (getline(infile, line)) {
        if (line.empty()) continue;

        SNPEntry e = parseLine(line);
        if (!e.valid) continue;

        const int informative_reads = e.ALT_ALT + e.ALT_REF + e.REF_ALT;
        if (informative_reads > 0) {
            coverages.push_back(informative_reads);
        }
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
            "HAP1=" + e.hap1 + " "
            "HAP2=" + e.hap2 + " "
            "HP_READS1=" + e.hp_reads1 + " "
            "HP_READS2=" + e.hp_reads2 + " "
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
            if (!divergentSameHp || hasSameAssignedHaplotype(e)) {
                f_div << formatted << '\n';
                retainedDivergent++;
            } else {
                filteredDivergent++;
            }
        } else if (cls.type == "snp1_before_snp2") {
            f_snp1 << formatted << '\n';
        } else if (cls.type == "snp2_before_snp1") {
            f_snp2 << formatted << '\n';
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
    f_err.close();

    if (divergentSameHp) {
        cout << "[divergent_same_hp] retained=" << retainedDivergent
             << " filtered=" << filteredDivergent << "\n";
    }
    cout << "Processing complete.\n";
}
