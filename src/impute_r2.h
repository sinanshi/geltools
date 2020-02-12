#ifndef IMPUTER2_H
#define IMPUTER2_H
#include "io.h"

class ImputeRsquare {

  public:
    ImputeRsquare() {}; 
    ImputeRsquare(Dosage& imp, Dosage& tru){
      impute_dose = imp;
      true_dose = tru;
    };
  Dosage impute_dose, true_dose;
  vector<int> allele_frequency_category;

  float r2(const vector<float>& imp, const vector<float>& tru);
  void read_allele_frequency_category(string afname);
  void aggregate_r2(string afname);

};


#endif
