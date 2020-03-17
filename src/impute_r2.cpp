#include "impute_r2.h"
#include <numeric>
#include <math.h>
#include <algorithm>

using namespace std;

float ImputeRsquare::r2(const vector<float>& est, const vector<float>& tru) {
    assert(est.size() == tru.size());
    int size = est.size();

/*    float diff;
    float truth_mean = accumulate(tru.begin(), tru.end(), 0.0) / size;
    float ss_total = 0.0;
    for (int i = 0; i < size; ++i) {
        diff = tru[i] - truth_mean;
        ss_total += diff * diff;
    };
    float ss_residual = 0.0;
    for (int i = 0; i < size; ++i) {
        diff = tru[i] - est[i];
        ss_residual += diff * diff;
    }
    return(1 - ss_residual / ss_total);*/

    float truth_mean = accumulate(tru.begin(), tru.end(), 0.0) / size;
    float est_mean = accumulate(est.begin(), est.end(), 0.0) / size;

    vector<float> diff_est(size);
    vector<float> diff_truth(size);

    for (int i = 0; i < size; ++i) diff_est[i] = est[i] - est_mean;
    for (int i = 0; i < size; ++i) diff_truth[i] = tru[i] - truth_mean;

    float cov = 0;
    float sd1 = 0;
    float sd2 = 0;
    for (int i = 0; i < size; ++i) {
        sd1 += diff_est[i] * diff_est[i];
    }; sd1 = sqrt(sd1 / size);
    for (int i = 0; i < size; ++i) {
        sd2 += diff_truth[i] * diff_truth[i];
    }; sd2 = sqrt(sd2 / size);

    for (int i = 0; i < size; ++i) {
        cov += diff_truth[i] * diff_est[i] ;
    };
//    cov /= size;
    cov = cov / sd1 / sd2 / size;

    return(cov * cov);

}

void ImputeRsquare::read_allele_frequency_category(string afname) {
    string line;
    ifstream affile(afname);
    if (affile.is_open()) {
        while(getline(affile, line)) {
            allele_frequency_category.push_back(stoi(line));
        }
        affile.close();
    } else
        throw(runtime_error(string("Cannot open the AF file ") + afname));
}

void ImputeRsquare::aggregate_r2(string afname) {
    read_allele_frequency_category(afname);
    assert((int)allele_frequency_category.size() == impute_dose.nsite());
    int n_category = *max_element(allele_frequency_category.begin(),
            allele_frequency_category.end());
    vector<vector<float>> impute_catg_dose(n_category + 1);
    vector<vector<float>> true_catg_dose(n_category + 1);

    for (int i = 0; i < impute_dose.nsite(); ++i) {
        if (allele_frequency_category[i] >= 0) {
            int k = allele_frequency_category[i];
            impute_catg_dose[k].insert(
                    impute_catg_dose[k].end(),
                    impute_dose.data[i].begin(), impute_dose.data[i].end());
            auto td = true_dose[impute_dose.rsid[i]];
            true_catg_dose[k].insert(true_catg_dose[k].end(),
                    td.begin(), td.end());
        }
    }

    ofstream r2file;
    r2file.open(afname + ".r2");
    for (int i = 0; i < n_category + 1; ++i) {
        r2file << r2(impute_catg_dose[i], true_catg_dose[i]) << endl;
    }
    r2file.close();
}
