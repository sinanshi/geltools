#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() 
#define CATCH_CONFIG_FAST_COMPILE
#include "catch.hpp"
#include "io.h"


TEST_CASE("Read VCF file - get the number of samples") {
  string fname = "../../tests/test.vcf.gz";
  VCFProc vcf_reader = VCFProc(fname);
  CHECK(vcf_reader.nsample == 3);
}


TEST_CASE("Read VCF file - read genotypes") {
  string fname = "../../tests/test.vcf.gz";
  VCFProc vcf_reader = VCFProc(fname);
  vcf_reader.read_all_genotype();

  CHECK(vcf_reader.geno.size() == 6);
  CHECK(vcf_reader.geno[0].size() == 6);


}

TEST_CASE("Read VCF file - msprime output") {
  string fname = "../../tests/simulated_data.vcf.gz";
  VCFProc vcf_reader = VCFProc(fname);
  vcf_reader.read_all_genotype();
  CHECK(vcf_reader.nsample == 25);
  CHECK(vcf_reader.nsite == 10740);
  CHECK(vcf_reader.polymorphic_pos[10739] == 2998870);
}


TEST_CASE("process row by row - dummy callback") {
  class TestCallBack {
    public:
      int num;
      TestCallBack(int num_):num(num_){}
      void inc(VCFRow& r) {
        cout << num << endl;
        num++;
      }
  };

  auto tcb = TestCallBack(0);
  function<void(VCFRow&)> inc = std::bind(&TestCallBack::inc, &tcb, _1);
  string fname = "../../tests/test.vcf.gz";
  VCFProc vproc(fname);
  vproc.process_row_by_row(inc);
}


TEST_CASE("read dosage from VCF files") {
  string fname = "../../tests/beagle_out_samll.vcf";
  VCFProc vcf_reader = VCFProc(fname);
  Dosage impute_dosage(false);
  vcf_reader.read_all_dosage_by_site(impute_dosage);
  auto idose1 = impute_dosage[string("NAT-1:47_A_TT")];
  CHECK(idose1[0] == 1.95f); 
  CHECK(idose1[1] == 1.11f);
  auto idose2 = impute_dosage[string("IND-1:143_A_AT")];
  CHECK(idose2[0] == 0.02f); 
  CHECK(idose2[1] == 0.5f);

  auto idose3 = impute_dosage[string("SNP-1:170_A_T")];
  CHECK(idose3[0] == 0.0f); 
  CHECK(idose3[1] == 1.0f);

  auto idose4 = impute_dosage[string("SNP-1:236_A_T")];
  CHECK(idose4[0] == 0.0f); 
  CHECK(idose4[1] == 0.0f);

  auto idose5 = impute_dosage[string("Not found")];
  CHECK(idose5.size() == 0); 


  Dosage true_dosage(true);
  vcf_reader.read_all_dosage_by_site(true_dosage);
  auto tdose1 = true_dosage[string("NAT-1:47_A_TT")];
  CHECK(tdose1[0] == 2.f); 
  CHECK(tdose1[1] == 1.f);
  auto tdose2 = true_dosage[string("IND-1:143_A_AT")];
  CHECK(tdose2[0] == 0.f); 
  CHECK(tdose2[1] == 0.f);

  auto tdose3 = true_dosage[string("SNP-1:170_A_T")];
  CHECK(tdose3[0] == 0.f); 
  CHECK(tdose3[1] == 1.f);

  auto tdose4 = true_dosage[string("SNP-1:236_A_T")];
  CHECK(tdose4[0] == 0.0f); 
  CHECK(tdose4[1] == 0.0f);

  auto tdose5 = true_dosage[string("Not found")];
  CHECK(tdose5.size() == 0); 

}

TEST_CASE("read pedigree file") {
  string fname = "../../tests/test_duotrio.vcf.gz";
  VCFProc vproc = VCFProc(fname);
  auto fam = Family(fname + ".ped", vproc);
  CHECK(fam.ntrio == 2);
  CHECK(fam.nduo == 2);
  CHECK(fam.ntrio_ped == 3);
  CHECK(fam.nduo_ped == 1);
}

TEST_CASE("switch error rate") {
  string fname = "../../tests/test_duotrio2.vcf";
  VCFProc vproc(fname);
  auto fam = Family(fname + ".ped", vproc);
  function<void(VCFRow&)> swe = bind(&Family::detect_switches, &fam, _1);
  vproc.process_row_by_row(swe);
  fam.write_summary();
}
