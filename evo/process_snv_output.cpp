#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <map>

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

float calculateReliabilityForRelationship(const SNPEntry& e, const string& relType, int minReads) {
    if (e.TOTAL == 0) return 0.0f;
    
    // Convert to percentages
    float obs_RR = (float)e.REF_REF / e.TOTAL;
    float obs_AR = (float)e.ALT_REF / e.TOTAL;
    float obs_RA = (float)e.REF_ALT / e.TOTAL;
    float obs_AA = (float)e.ALT_ALT / e.TOTAL;
    
    // Expected maximum proportions based on VAFs
    float exp_RR = min(1.0f - e.vaf1, 1.0f - e.vaf2);
    float exp_AR = min(e.vaf1, 1.0f - e.vaf2);
    float exp_RA = min(1.0f - e.vaf1, e.vaf2);
    float exp_AA = min(e.vaf1, e.vaf2);
    
    float reliability = 0.0f;
    
    if (relType == "snp1_before_snp2") {
        // Required: R/R, A/R, A/A (each must have >= minReads)
        if (e.REF_REF < minReads || e.ALT_REF < minReads || e.ALT_ALT < minReads) {
            return 0.0f;
        }
        
        float min_ratio = 1.0f;
        if (exp_RR > 0) min_ratio = min(min_ratio, obs_RR / exp_RR);
        if (exp_AR > 0) min_ratio = min(min_ratio, obs_AR / exp_AR);
        if (exp_AA > 0) min_ratio = min(min_ratio, obs_AA / exp_AA);
        
        // Forbidden: R/A
        float forbidden_penalty = obs_RA;
        
        reliability = min_ratio - forbidden_penalty;
        
    } else if (relType == "snp2_before_snp1") {
        // Required: R/R, R/A, A/A
        if (e.REF_REF < minReads || e.REF_ALT < minReads || e.ALT_ALT < minReads) {
            return 0.0f;
        }
        
        float min_ratio = 1.0f;
        if (exp_RR > 0) min_ratio = min(min_ratio, obs_RR / exp_RR);
        if (exp_RA > 0) min_ratio = min(min_ratio, obs_RA / exp_RA);
        if (exp_AA > 0) min_ratio = min(min_ratio, obs_AA / exp_AA);
        
        // Forbidden: A/R
        float forbidden_penalty = obs_AR;
        
        reliability = min_ratio - forbidden_penalty;
        
    } else if (relType == "cooccurring") {
        // Required: R/R, A/A
        if (e.REF_REF < minReads || e.ALT_ALT < minReads) {
            return 0.0f;
        }
        
        float min_ratio = 1.0f;
        if (exp_RR > 0) min_ratio = min(min_ratio, obs_RR / exp_RR);
        if (exp_AA > 0) min_ratio = min(min_ratio, obs_AA / exp_AA);
        
        // Forbidden: A/R, R/A
        float forbidden_penalty = obs_AR + obs_RA;
        
        reliability = min_ratio - forbidden_penalty;
        
    } else if (relType == "cooccurring_loss") {
        // Required: A/A only
        if (e.ALT_ALT < minReads) {
            return 0.0f;
        }
        
        float min_ratio = 1.0f;
        if (exp_AA > 0) min_ratio = min(min_ratio, obs_AA / exp_AA);
        
        // Forbidden: R/R, A/R, R/A
        float forbidden_penalty = obs_RR + obs_AR + obs_RA;
        
        reliability = min_ratio - forbidden_penalty;
        
    } else if (relType == "divergent") {
        // Required: A/R, R/A
        if (e.ALT_REF < minReads || e.REF_ALT < minReads) {
            return 0.0f;
        }
        
        float min_ratio = 1.0f;
        if (exp_AR > 0) min_ratio = min(min_ratio, obs_AR / exp_AR);
        if (exp_RA > 0) min_ratio = min(min_ratio, obs_RA / exp_RA);
        
        // Forbidden: R/R, A/A
        float forbidden_penalty = obs_RR + obs_AA;
        
        reliability = min_ratio - forbidden_penalty;
        
    } else if (relType == "snp1_before_snp2_loss") {
        // Required: A/R, A/A
        if (e.ALT_REF < minReads || e.ALT_ALT < minReads) {
            return 0.0f;
        }
        
        float min_ratio = 1.0f;
        if (exp_AR > 0) min_ratio = min(min_ratio, obs_AR / exp_AR);
        if (exp_AA > 0) min_ratio = min(min_ratio, obs_AA / exp_AA);
        
        // Forbidden: R/R, R/A
        float forbidden_penalty = obs_RR + obs_RA;
        
        reliability = min_ratio - forbidden_penalty;
        
    } else if (relType == "snp2_before_snp1_loss") {
        // Required: R/A, A/A
        if (e.REF_ALT < minReads || e.ALT_ALT < minReads) {
            return 0.0f;
        }
        
        float min_ratio = 1.0f;
        if (exp_RA > 0) min_ratio = min(min_ratio, obs_RA / exp_RA);
        if (exp_AA > 0) min_ratio = min(min_ratio, obs_AA / exp_AA);
        
        // Forbidden: R/R, A/R
        float forbidden_penalty = obs_RR + obs_AR;
        
        reliability = min_ratio - forbidden_penalty;
    }
    
    return reliability;
}

// Classify based on reliability scores
struct Classification {
    string type;
    float final_reliability;
    float best_score;
    float margin;
};

Classification classifyEntry(const SNPEntry& e, int minReads) {
    // Calculate reliability for all relationship types
    map<string, float> reliabilities;
    reliabilities["snp1_before_snp2"] = calculateReliabilityForRelationship(e, "snp1_before_snp2", minReads);
    reliabilities["snp2_before_snp1"] = calculateReliabilityForRelationship(e, "snp2_before_snp1", minReads);
    reliabilities["cooccurring"] = calculateReliabilityForRelationship(e, "cooccurring", minReads);
    reliabilities["cooccurring_loss"] = calculateReliabilityForRelationship(e, "cooccurring_loss", minReads);
    reliabilities["divergent"] = calculateReliabilityForRelationship(e, "divergent", minReads);
    reliabilities["snp1_before_snp2_loss"] = calculateReliabilityForRelationship(e, "snp1_before_snp2_loss", minReads);
    reliabilities["snp2_before_snp1_loss"] = calculateReliabilityForRelationship(e, "snp2_before_snp1_loss", minReads);
    
    // Find best and second-best
    string best_type = "error";
    float best_score = -999.0f;
    float second_best = -999.0f;
    
    for (const auto& [type, score] : reliabilities) {
        if (score > best_score) {
            second_best = best_score;
            best_score = score;
            best_type = type;
        } else if (score > second_best) {
            second_best = score;
        }
    }
    
    // Calculate margin and final reliability
    float margin = best_score - second_best;
    float final_reliability = best_score * (1.0f + margin);
    
    // If best score is <= 0, mark as error
    if (best_score <= 0.0f) {
        best_type = "error";
        final_reliability = 0.0f;
    }
    
    Classification result;
    result.type = best_type;
    result.final_reliability = final_reliability;
    result.best_score = best_score;
    result.margin = margin;
    
    return result;
}

void processFile(const string& filename, int minReads) {
    ifstream infile(filename);
    if (!infile) {
        cerr << "Error: cannot open input file.\n";
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

    string line;
    while (getline(infile, line)) {
        if (line.empty()) continue;

        SNPEntry e = parseLine(line);
        if (!e.valid) continue;

        Classification cls = classifyEntry(e, minReads);
        
        // Format output with reliability scores
        string formatted =
            e.chr + " " + to_string(e.pos1) + " " + to_string(e.pos2) + " "
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