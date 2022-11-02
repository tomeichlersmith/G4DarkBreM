/** 
 * @file extract-library.cxx
 * definition of g4db-extract-library executable
 */

#include <iostream>
#include <fstream>

#include "G4DarkBreM/ParseLibrary.h"

/**
 * printout how to use g4db-extract-library
 */
void usage() {
  std::cout << 
      "USAGE:\n"
      "  g4db-extract-library [options] db-lib\n"
      "\n"
      "  Extract the input DB event library into a single CSV file\n"
      "\n"
      "ARGUMENTS\n"
      "  db-lib : dark brem event library to load and extract\n"
      "\n"
      "OPTIONS\n"
      "  -h,--help             : produce this help and exit\n"
      "  -o,--output           : output file to write extracted events to\n"
      "                          use the input library name with the '.csv' extension added by default\n"
      "  --aprime-id           : A' ID number as used in the LHE files\n"
      << std::flush;
}

/**
 * definition of g4db-extract-library
 */
int main(int argc, char* argv[]) try {
  std::string db_lib{};
  std::string output_filename{};
  int aprime_id{622};
  for (int i_arg{1}; i_arg < argc; i_arg++) {
    std::string arg{argv[i_arg]};
    if (arg == "-h" or arg == "--help") {
      usage();
      return 0;
    } else if (arg == "-o" or arg == "--output") {
      if (i_arg+1 >= argc) {
        std::cerr << arg << " requires an argument after it" << std::endl;
        return 1;
      }
      output_filename = argv[++i_arg];
    } else if (arg == "--aprime-id") {
      if (i_arg+1 >= argc) {
        std::cerr << arg << " requires an argument after it" << std::endl;
        return 1;
      }
      aprime_id = std::stoi(argv[++i_arg]);
    } else if (not arg.empty() and arg[0] == '-') {
      std::cerr << arg << " is not a recognized option" << std::endl;
      return 1;
    } else {
      db_lib = arg;
    }
  }

  if (db_lib.empty()) {
    std::cerr << "ERROR: DB event library not provided." << std::endl;
    return 1;
  }

  if (output_filename.empty()) {
    // remove trailing slash if present
    if (db_lib.back() == '/') db_lib.pop_back();
    output_filename = db_lib+".csv";
  }

  std::ofstream output{output_filename};
  if (not output.is_open()) {
    std::cerr << "ERROR: Unable to open " << output_filename << " for writing." << std::endl;
    return 2;
  }

  std::map<double, std::vector<g4db::OutgoingKinematics>> lib;
  parseLibrary(db_lib, aprime_id, lib);
  dumpLibrary(output, lib);

  output.close();
  return 0;
} catch (const std::exception& e) {
  std::cerr << "ERROR: " << e.what() << std::endl;
  return 127;
}
