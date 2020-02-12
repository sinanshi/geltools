#include "parameters.h"

Options::Options(int argc, char **argv) {
  bpo::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
    ("mode", bpo::value<string>()->required(), "mode switch")
    ("input", bpo::value<string>()->required(), "vcf input file path")
    ("truth", bpo::value<string>()->required(), "vcf truth file path")
    ("freq", bpo::value<string>()->required(), "Allele frequency file path")
    ("no_dosage", bpo::value<string>(), "Allele frequency file path");

  bpo::variables_map vm;
  bpo::store(bpo::parse_command_line(argc, argv, desc), vm);
  bpo::notify(vm);
  
  if (vm.count("help")) {
    cout << "help here" <<endl;
    exit(0);
  }

  if (vm.count("mode"))
      mode = vm["mode"].as<string>();
  if (vm.count("input"))
      input_vcf_file_path = vm["input"].as<string>();

  if (mode == "r2") {
      if (vm.count("truth")) truth_file_path = vm["truth"].as<string>();
      else throw(runtime_error("argument --truth is not specified."));

      if (vm.count("freq")) af_file_path = vm["freq"].as<string>();
      else throw(runtime_error("argument --freq is not specified."));
      
      if (vm.count("no_dosage")) no_dosage = true;
      else no_dosage = false;
  }

}
