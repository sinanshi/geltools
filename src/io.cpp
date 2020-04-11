#include "io.h"
#include <cstring>
#include <stdexcept>
#include <assert.h>
#include <boost/tokenizer.hpp>


typedef boost::tokenizer< boost::escaped_list_separator<char> > Tokenizer;

htsFile* VCFProc::open_bcf_file(string fn){
  char fchar[fn.size() + 1];
  strcpy(fchar, fn.c_str());
  htsFile *file = bcf_open(fchar, "r");
  if (file == NULL)
    throw(runtime_error(fn + " cannot be found!\n"));
  return(file);
}


VCFProc::VCFProc(string file_path) {
  fname = file_path;
  auto file = open_bcf_file(file_path);
  bcf_hdr_t *hdr = bcf_hdr_read(file);
  nsample = bcf_hdr_nsamples(hdr);
  nseq = nsample * 2;
  for (int i=0; i < nsample; ++i) {
   samples_map.insert(pair<string, Sample>(hdr->samples[i], Sample(hdr->samples[i], i)));
  }
  for (int i = 0; i < nsample * 2; ++i) {
    geno.push_back(vector<bool>());
  }
  bcf_hdr_destroy(hdr);
  bcf_close(file);
}


void VCFProc::read_all_genotype() {
  nsite = 0;
  int32_t *gt =NULL;
  int ngt_arr = 0;

  polymorphic_pos.clear();

  htsFile *file = open_bcf_file(fname);
  bcf_hdr_t *hdr = bcf_hdr_read(file);
  bcf1_t *rec = bcf_init();
  while(bcf_read(file, hdr, rec) == 0) {
    bcf_get_genotypes(hdr, rec, &gt, &ngt_arr);
    assert(ngt_arr == nsample * 2);
    polymorphic_pos.push_back(rec->pos + 1); // don't know why the pos from htslib is always 1 number off from the reading.
    for (int i = 0; i < nsample * 2; ++i) {
      //! TODO: here treat all the missing or non-biallelic
      //! case as 0. This should be changed.
      geno[i].push_back(bcf_gt_allele(gt[i]) == 1);
    }
  }
  if (geno[0].size() == 0)
    throw(runtime_error("VCFreader doesn't contain any sample or the reader was not initialized properly.\n"));
  
  nsite = polymorphic_pos.size();
  free(gt);
  bcf_close(file);
  bcf_destroy(rec);
}

void VCFProc::read_all_dosage_by_site(Dosage& dosage) {
  nsite = 0;
  float *ds = NULL;
  int nds_arr = 0;

  int32_t *gt = NULL;
  int ngt_arr = 0;

  int32_t *imp = NULL;
  int nimp = 0;

  polymorphic_pos.clear();

  htsFile *file = open_bcf_file(fname);
  bcf_hdr_t *hdr = bcf_hdr_read(file);
  bcf1_t *rec = bcf_init();
  string variant_type;
  dosage.nimputed = 0;
  cout << "reading dosages (" << nsample << ")..." << endl;
  while(bcf_read(file, hdr, rec) == 0) {
    bcf_unpack(rec, BCF_UN_STR);
    //! remove all the none bi-allelic cases.
    if (rec-> n_allele == 2) {
      polymorphic_pos.push_back(rec->pos + 1);
      bool is_imputed = bcf_get_info_flag(hdr, rec, "IMP", &imp, &nimp);

      if (is_imputed) {
        dosage.site_is_imputed.push_back(true);
        dosage.nimputed++;
      } else 
        dosage.site_is_imputed.push_back(false);

      if (dosage.genotype_as_dosage) {
        bcf_get_genotypes(hdr,rec, &gt, &ngt_arr);
        assert(ngt_arr == nsample * 2);
        dosage.data.push_back(vector<float>(nsample));
        for (int i = 0; i < nsample; ++i) {
          dosage.data.back()[i] = bcf_gt_allele(gt[2 * i]) + bcf_gt_allele(gt[2 * i + 1]);
        }
      }else {
        bcf_get_format_float(hdr, rec, "DS", &ds, &nds_arr);
        assert(nds_arr == nsample);
        dosage.data.push_back(vector<float>(ds, ds + nsample));
      }

      // id is in the form of TYPE-chr:pos_ref_alt
      string id = get_vcf_row_id(rec, hdr);
      dosage.id_row_map[id] = dosage.data.size() - 1;
      dosage.rsid.push_back(id);
    }
      nmultiallelic++;
  }
  nsite = polymorphic_pos.size();
  free(ds);
  bcf_close(file);
  bcf_destroy(rec);
}

string VCFProc::get_vcf_row_id(bcf1_t* rec, bcf_hdr_t *hdr) {
  string vartype = "NAT";
  int vcftype = bcf_get_variant_types(rec);
  if (vcftype == VCF_SNP) vartype = "SNP";
  if (vcftype == VCF_INDEL) vartype = "IND";

  string id =  vartype + string("-") + string(bcf_hdr_id2name(hdr, rec->rid)) + 
    string(":") + to_string(rec->pos + 1) + string("_") + string(rec->d.allele[0]) + 
    string("_") + string(rec->d.allele[1]);
  return(id);
}


void Family::read_pedigree(const VCFProc& vcf) {
  cout << "processing the pedigree ...";
  int nfather, nmother, nchild;
  Sample child, father, mother;
  string line;
  vector<string> vec;
  ifstream data_in(ped_path.c_str());
  if (!data_in.is_open()) throw(runtime_error("The pedigree file cannot be found.\n"));

  while(getline(data_in, line)) {
    nchild = nmother = nfather = 0;
    Tokenizer tok(line);
    vec.assign(tok.begin(), tok.end());
    if(vec.size() != 3)
      throw(runtime_error("Pedigree file doesn't contain 3 columns\n"));

    if (vec[1] != "0" && vec[2] != "0") ntrio_ped++;
    if (vec[1] == "0" || vec[2] == "0") nduo_ped++;
    if (vec[1] == "0" && vec[2] == "0")// we can include this case later.
      throw(runtime_error("pedigree cannot contain unrealted sample.\n"));
    if (vcf.samples_map.find(vec[0]) != vcf.samples_map.end()) {
      child = vcf.samples_map.at(vec[0]); nchild++;
    }
    if (vcf.samples_map.find(vec[1]) != vcf.samples_map.end()) {
      father = vcf.samples_map.at(vec[1]); nfather++;
    }
    if (vcf.samples_map.find(vec[2]) != vcf.samples_map.end()) {
      mother = vcf.samples_map.at(vec[2]); nmother++;
    }

    if (nfather + nmother + nchild == 0) nmiss++;
    else {
      if (nfather + nmother == 2) { // trios
        trios.push_back(Trio(child, father, mother));

      } else { // duos
        if (nfather == 1)  // duo with father and child
          duos.push_back(Duo(child, father));
        else  // duo with mother and child
          duos.push_back(Duo(child, mother));
        }
      }
    }
  ntrio = trios.size();
  nduo = duos.size();
  cout << endl << ntrio_ped << " trios and " << nduo_ped << 
    " duos have been found in the pedigree file." << endl;
  cout << ntrio << " trios and " << nduo << " duos " << 
    "have been detected from the VCF file" << endl;
}

void Trio::update_switch_state(VCFRow& vr) {
  int  c0 = vr.get_sample_left_gt(child); // child left strand
  int  c1 = vr.get_sample_right_gt(child); // right strand
  bool child_is_phased = vr.get_sample_phase(child); // child is phased?
  int  f0 = vr.get_sample_left_gt(father);//father left strand
  int  f1 = vr.get_sample_right_gt(father);//father right strand
  bool father_is_phased = vr.get_sample_phase(father);//father is phased?
  int  m0 = vr.get_sample_left_gt(mother);//mother left strand
  int  m1 = vr.get_sample_right_gt(mother);//mother right strand
  bool mother_is_phased = vr.get_sample_phase(mother);//is mother phased?

  int cg = c0 + c1;
  int fg = f0 + f1;
  int mg = m0 + m1;

  bool cphase, fphase, mphase, ctest, ftest, mtest;
  ctest = ftest = mtest = 0;

  if ((fg == 0 && mg == 0 && cg == 1) || (fg == 0 && mg == 0 && cg == 2) ||
      (fg == 0 && mg == 1 && cg == 2) || (fg == 0 && mg == 2 && cg == 0) ||
      (fg == 0 && mg == 2 && cg == 2) || (fg == 1 && mg == 2 && cg == 0) ||
      (fg == 1 && mg == 0 && cg == 2) || (fg == 2 && mg == 0 && cg == 2) ||
      (fg == 2 && mg == 1 && cg == 0) || (fg == 2 && mg == 2 && cg == 0) ||
      (fg == 2 && mg == 2 && cg == 1) || (fg == 2 && mg == 0 && cg == 0)) {
    child.nmendel++; father.nmendel++; mother.nmendel++;
  } 
  else if (fg == 0 && mg == 1 && cg == 0) { fphase = (f0 == 0 && f1 == 0); mphase = (m0 == 1 && m1 == 0); cphase = (c0 == 0 && c1 == 0); mtest=1;}
  else if (fg == 0 && mg == 1 && cg == 1) { fphase = (f0 == 0 && f1 == 0); mphase = (m0 == 0 && m1 == 1); cphase = (c0 == 0 && c1 == 1); mtest=1; ctest=1;}
  else if (fg == 0 && mg == 2 && cg == 1) { fphase = (f0 == 0 && f1 == 0); mphase = (m0 == 1 && m1 == 1); cphase = (c0 == 0 && c1 == 1); ctest=1;}
  else if (fg == 1 && mg == 0 && cg == 0) { fphase = (f0 == 0 && f1 == 1); mphase = (m0 == 0 && m1 == 0); cphase = (c0 == 0 && c1 == 0); ftest=1;}
  else if (fg == 1 && mg == 0 && cg == 1) { fphase = (f0 == 1 && f1 == 0); mphase = (m0 == 0 && m1 == 0); cphase = (c0 == 1 && c1 == 0); ftest=1; ctest=1;}
  else if (fg == 1 && mg == 1 && cg == 0) { fphase = (f0 == 0 && f1 == 1); mphase = (m0 == 1 && m1 == 0); cphase = (c0 == 0 && c1 == 0); ftest=1; mtest=1;}
  else if (fg == 1 && mg == 1 && cg == 2) { fphase = (f0 == 1 && f1 == 0); mphase = (m0 == 0 && m1 == 1); cphase = (c0 == 1 && c1 == 1); ftest=1; mtest=1;}
  else if (fg == 1 && mg == 2 && cg == 1) { fphase = (f0 == 0 && f1 == 1); mphase = (m0 == 1 && m1 == 1); cphase = (c0 == 0 && c1 == 1); ftest=1; ctest=1;}
  else if (fg == 1 && mg == 2 && cg == 2) { fphase = (f0 == 1 && f1 == 0); mphase = (m0 == 1 && m1 == 1); cphase = (c0 == 1 && c1 == 1); ftest=1; }
  else if (fg == 2 && mg == 0 && cg == 1) { fphase = (f0 == 1 && f1 == 1); mphase = (m0 == 0 && m1 == 0); cphase = (c0 == 1 && c1 == 0); ftest=1; ctest=1;}
  else if (fg == 2 && mg == 1 && cg == 1) { fphase = (f0 == 1 && f1 == 1); mphase = (m0 == 1 && m1 == 0); cphase = (c0 == 1 && c1 == 0); mtest=1; ctest=1;}
  else if (fg == 2 && mg == 1 && cg == 2) { fphase = (f0 == 1 && f1 == 1); mphase = (m0 == 0 && m1 == 1); cphase = (c0 == 1 && c1 == 1); mtest=1;
  }


  if (ctest && child_is_phased) {
    if (cphase != child.current_phase) {
      child.nswitch++;
      child.ntest++;
      child.current_phase = !child.current_phase;
    }
  }

  if (ftest && father_is_phased) {
    if (fphase != father.current_phase) {
      father.nswitch++;
      father.ntest++;
      father.current_phase = !father.current_phase;
    }
  }

  if (mtest && mother_is_phased) {
    if (mphase != mother.current_phase) {
      mother.nswitch++;
      mother.ntest++;
      mother.current_phase = !mother.current_phase;
    }
  }
}


void Family::detect_switches(VCFRow & vr) {
  for (size_t t = 0; t < trios.size(); ++t)
    trios[t].update_switch_state(vr);
}

void Family::write_summary() {
  ofstream swefile;
  swefile.open(ped_path + ".swe", ios::out);
  swefile << "type\tchild\tfather\tmother\tntest\tnmendel\tnswitch" << endl;
  for (auto t : trios) {
    swefile << "trios\t" << t.child.name << "\t" << t.father.name << "\t" <<
      t.mother.name << "\t" << t.child.ntest << "\t" << t.child.nmendel << "\t" <<
      t.child.nswitch << endl;
  }
  swefile.close();
}
