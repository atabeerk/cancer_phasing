#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

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

    // Store original (pre-threshold) counts
    int orig_ALT_ALT = 0;
    int orig_ALT_REF = 0;
    int orig_REF_ALT = 0;
    int orig_REF_REF = 0;
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
        int val = stoi(kv.substr(pos + 1));
        if (key == "ALT_ALT") e.ALT_ALT = val;
        else if (key == "ALT_REF") e.ALT_REF = val;
        else if (key == "REF_ALT") e.REF_ALT = val;
        else if (key == "REF_REF") e.REF_REF = val;
        else if (key == "TOTAL") e.TOTAL = val;
    }

    // TOTAL must exist and be > 0
    if (e.TOTAL <= 0) e.valid = false;
    return e;
}

void applyThreshold(SNPEntry& e) {
    // Store originals first
    e.orig_ALT_ALT = e.ALT_ALT;
    e.orig_ALT_REF = e.ALT_REF;
    e.orig_REF_ALT = e.REF_ALT;
    e.orig_REF_REF = e.REF_REF;

    double threshold = max((e.TOTAL * 0.12), 2.0);
    if (e.ALT_ALT < threshold) e.ALT_ALT = 0;
    if (e.ALT_REF < threshold) e.ALT_REF = 0;
    if (e.REF_ALT < threshold) e.REF_ALT = 0;
    if (e.REF_REF < threshold) e.REF_REF = 0;
}

void processFile(const string& filename) {
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
    ofstream f_err(base + "_errors.txt");

    string line;
    while (getline(infile, line)) {
        if (line.empty()) continue;

        SNPEntry e = parseLine(line);
        if (!e.valid) continue;  // skip TOTAL <= 0

        applyThreshold(e);

        // --- Create formatted output including ORIGINAL counts ---
        string formatted =
            e.chr + " " + to_string(e.pos1) + " " + to_string(e.pos2) + " "
            "ALT_ALT=" + to_string(e.ALT_ALT) + " "
            "ALT_REF=" + to_string(e.ALT_REF) + " "
            "REF_ALT=" + to_string(e.REF_ALT) + " "
            "REF_REF=" + to_string(e.REF_REF) + " "
            "TOTAL=" + to_string(e.TOTAL) + " "
            "ORIGINAL=" + to_string(e.orig_ALT_ALT) + "/" + to_string(e.orig_ALT_REF) + "/" +
                          to_string(e.orig_REF_ALT) + "/" + to_string(e.orig_REF_REF);

        bool alt_alt = e.ALT_ALT > 0;
        bool alt_ref = e.ALT_REF > 0;
        bool ref_alt = e.REF_ALT > 0;
        bool ref_ref = e.REF_REF > 0;

        if (alt_alt && ref_ref && !alt_ref && !ref_alt)
            f_co << formatted << '\n';                                // Co-occurring
        else if (!alt_alt && alt_ref && ref_alt)
            f_div << formatted << '\n';                               // Divergent
        else if (alt_ref && alt_alt && ref_ref && !ref_alt)
            f_snp1 << formatted << '\n';                              // SNP1 before SNP2
        else if (ref_alt && alt_alt && ref_ref && !alt_ref)
            f_snp2 << formatted << '\n';                              // SNP2 before SNP1
        else if (alt_alt && alt_ref && !ref_alt && !ref_ref)
            f_snp1_loss << formatted << '\n';                         // SNP1 before SNP2 (possible loss)
        else if (alt_alt && ref_alt && !alt_ref && !ref_ref)
            f_snp2_loss << formatted << '\n';                         // SNP2 before SNP1 (possible loss)
        else
            f_err << formatted << '\n';                               // Ambiguous
    }

    infile.close();
    f_co.close();
    f_div.close();
    f_snp1.close();
    f_snp2.close();
    f_snp1_loss.close();
    f_snp2_loss.close();
    f_err.close();

    cout << "Processing complete.\n";
}

