#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() 
#define CATCH_CONFIG_FAST_COMPILE
#include "catch.hpp"
#include "impute_r2.h"


TEST_CASE("Compare the true value of R^2") {
  auto ir2 = ImputeRsquare();
  vector<float> x1 = {1,2,4};
  vector<float> y1 = {1,2,4};

  CHECK(ir2.r2(x1, y1) == 1.f);

  vector<float> y = {0.841471,0.9092974,0.14112,-0.7568025,
    -0.9589243,-0.2794155,0.6569866,0.9893582,0.4121185,-0.5440211};
  vector<float> yhat = {0.32623,0.28509,0.24395,0.20281,0.16167,
    0.12053,0.07939,0.03825,-0.00289,-0.04403};
  CHECK(ir2.r2(yhat, y) == Approx(0.02907).epsilon(0.0001));
}

TEST_CASE("Calculate the aggregate R2") {

  string fname = "../../tests/beagle_out2.vcf";
  string afname = "../../tests/beagle_out2.af";
  VCFProc vcf_reader = VCFProc(fname);
  Dosage impute_dosage(false);
  Dosage true_dosage(true);
  vcf_reader.read_all_dosage_by_site(impute_dosage);
  vcf_reader.read_all_dosage_by_site(true_dosage);

  auto rsq = ImputeRsquare(impute_dosage, true_dosage);
  rsq.aggregate_r2(afname);

}
