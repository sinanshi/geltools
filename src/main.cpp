#include <iostream>
#include <Eigen/Dense>
#include "parameters.h"
#include "io.h"
#include "impute_r2.h"

int main(int argc, char **argv)
{

  Options params(argc, argv);
  if (params.mode == "switch") {
    VCFProc vcf(params.input_vcf_file_path);
    auto fam = Family(params.input_vcf_file_path + ".ped", vcf);
    function<void(VCFRow&)> swe = bind(&Family::detect_switches, &fam, _1);
    vcf.process_row_by_row(swe);
    fam.write_summary();
  }

  if (params.mode == "r2") {
      VCFProc vcf_reader_impute = VCFProc(params.input_vcf_file_path);
      Dosage impute_dosage;
      if (params.no_dosage)
          impute_dosage.genotype_as_dosage = true;
      else
          impute_dosage.genotype_as_dosage = false;

      VCFProc vcf_reader_truth = VCFProc(params.truth_file_path);
      Dosage true_dosage(true);

      vcf_reader_impute.read_all_dosage_by_site(impute_dosage);
      vcf_reader_truth.read_all_dosage_by_site(true_dosage);

      ImputeRsquare rsq(impute_dosage, true_dosage);
      rsq.aggregate_r2(params.af_file_path);
  }
  return 0;
}

