#ifndef PARAMETERS_H
#define PARAMETERS_H
#include <iostream>
#include <string>
#include <boost/program_options.hpp>

using namespace std;
namespace bpo = boost::program_options;

//! This is for the parameters of EM algorithm.
class Options {
  public:
  string input_vcf_file_path;
  string af_file_path;
  string truth_file_path;
  bool no_dosage;
  string mode;
    Options(int argc, char **argv);
};

#endif
