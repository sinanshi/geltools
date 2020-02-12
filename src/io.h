#ifndef IO_H
#define IO_H
#include "htslib/vcf.h"
#include <vector>
#include <string>
#include <iostream>
#include <memory>
#include <fstream>
#include <functional>
#include <map>
using namespace std; 
using namespace std::placeholders;


class Sample {
  public:
    Sample(){};
    int nswitch = 0;
    int nflip = 0;
    int nmendel = 0;
    int ntest = 0;
    int nA, nG, nC, nT;
    string name;
    int ind;
    bool current_phase = true;
    //! @param sample name
    //! @param index within the VCF file
    Sample(string nm, int index): name(nm), ind(index){};
};

class Dosage {
    public:
        vector<float> & operator[](string id) { 
            std::map<string, int>::iterator id_iter = id_row_map.find(id);
            if (id_iter != id_row_map.end()) 
                return(data[id_iter->second]);
            else {
                return nulldose;
            }
        }
        Dosage(bool gasd) {genotype_as_dosage = gasd;}
        Dosage(){};
        int nsite() {return(data.size()); }
        int nsample() {return(data[0].size()); }
        bool genotype_as_dosage;
        map<string, int> id_row_map;
        vector<string> rsid;
        vector<vector<float>> data;
    private:
        vector<float> nulldose = vector<float>();
};



class VCFRow {
  public:
    int32_t *gt = NULL;
    int get_sample_left_gt(Sample& s) {
      int32_t hap = bcf_gt_allele(gt[s.ind * 2]);
      if (gt[s.ind * 2] == bcf_gt_missing) return(-2);
      if (gt[s.ind * 2] == bcf_int32_vector_end) 
        throw(runtime_error("something is wrong while parsing the genotype from VCF.\n"));
      return(hap);
    };

    int get_sample_right_gt(Sample& s) {
      int32_t hap = bcf_gt_allele(gt[s.ind * 2 + 1]);
      if (gt[s.ind * 2 + 1] == bcf_gt_missing) return(-2);
      if (gt[s.ind * 2 + 1] == bcf_int32_vector_end) 
        throw(runtime_error("something is wrong while parsing the genotype from VCF.\n"));
      return(hap);
    };

    bool get_sample_phase(Sample &s) {
      return(bcf_gt_is_phased(gt[s.ind * 2 + 1]));
    }

    VCFRow(){};
    ~VCFRow() {
      free(gt);
    }
};

class VCFProc {
  public:
    string fname;
    map<string, Sample> samples_map;
    int nsample, nseq, nsite, nmultiallelic;
    VCFProc(string file_path);
    //! The sites recorded in the VCF file.
    vector<float> polymorphic_pos;

    vector<vector<bool>> geno;
    void read_all_genotype();
    void read_all_dosage_by_site(Dosage& dosage);
    string get_vcf_row_id(bcf1_t* rec, bcf_hdr_t *hdr);
    void process_row_by_row(function<void(VCFRow&)> &func) {
      VCFRow row;
      int ngt_arr = 0;
      htsFile *file = open_bcf_file(fname);
      bcf_hdr_t *hdr = bcf_hdr_read(file);
      bcf1_t *rec = bcf_init();
      while(bcf_read(file, hdr, rec) == 0) {
        bcf_get_genotypes(hdr, rec, &row.gt, &ngt_arr);
        func(row);
      }
      bcf_close(file);
      bcf_destroy(rec);
    };
  private:
    //! This is just a wrapper of bcf_open in htslib
    //! but taking string as an input.
    htsFile *open_bcf_file(string fn);

};

class Trio {
  public:
    Trio(Sample c, Sample f, Sample m) : child(c), father(f), mother(m) {};
    Sample child, father, mother;
    void update_switch_state(VCFRow &vr);
};

class Duo {
  public:
    Duo(Sample c, Sample p): child(c), parent(p) {};
    Sample child, parent;
};




class Family {
  public:
    string ped_path;
    Family(string pp,const VCFProc& vcf):ped_path(pp) {
      ntrio = nduo = nmiss = ntrio_ped = nduo_ped = 0;
      read_pedigree(vcf);
    };
    int ntrio, nduo, nmiss;
    int ntrio_ped, nduo_ped;
    vector<Trio> trios;
    vector<Duo> duos;
    //! Detect switches, flips and mendel error.
    //! This function process every row of the bcf file, 
    //! hence need to be fed to VCFProc.process_row_by_row.
    //! Note: we can only deal with the bi-allelic case!
    void detect_switches(VCFRow&);
    void write_summary();
  private:
    void read_pedigree(const VCFProc&);
};


#endif
