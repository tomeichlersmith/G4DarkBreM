/**
 * @file xsec_calc.cxx
 * definition of g4db-xsec-calc executable
 */

#include <fstream>
#include <iostream>
#include <unistd.h>

#include "G4DarkBreM/G4DarkBremsstrahlung.h"
#include "G4DarkBreM/G4DarkBreMModel.h"
#include "G4DarkBreM/G4APrime.h"

/**
 * print out how to use g4db-xsec-calc
 */
void usage() {
  std::cout <<
    "USAGE:\n"
    "  g4db-xsec-calc [options]\n"
    "\n"
    "Calculate dark brem cross sections and write them out to a CSV table\n"
    "\n"
    "OPTIONS\n"
    "  -h,--help    : produce this help and exit\n"
    "  -o,--output  : output file to write scaled events to\n"
    "  -M,--ap-mass : mass of dark photon in GeV\n"
    "  --muons      : pass to set lepton to muons (otherwise electrons)\n"
    "  --energy     : python-like arange for input energies in GeV (stop, start stop, start stop step)\n"
    "                 default start is 0 and default step is 0.1 GeV\n"
    "  --target     : define target material with two parameters (atomic units): Z A\n"
    << std::flush;
}

/**
 * definition of g4db-xsec-calc
 *
 * We use the cross section caching table used within the G4DarkBremsstrahlung process.
 */
int main(int argc, char* argv[]) try {
  std::string output_filename{"xsec.csv"};
  double ap_mass{0.1};
  double min_energy{0.};
  double max_energy{4.};
  double energy_step{0.1};
  double target_Z{74.};
  double target_A{183.84};
  bool muons{false};
  for (int i_arg{1}; i_arg < argc; ++i_arg) {
    std::string arg{argv[i_arg]};
    if (arg == "-h" or arg == "--help") {
      usage();
      return 0;
    } else if (arg == "--muons") {
      muons = true;
    } else if (arg == "-o" or arg == "--output") {
      if (i_arg+1 >= argc) {
        std::cerr << arg << " requires an argument after it" << std::endl;
        return 1;
      }
      output_filename = argv[++i_arg];
    } else if (arg == "-M" or arg == "--ap-mass") {
      if (i_arg+1 >= argc) {
        std::cerr << arg << " requires an argument after it" << std::endl;
        return 1;
      }
      ap_mass = std::stod(argv[++i_arg]);
    } else if (arg == "--energy") {
      std::vector<std::string> args;
      while (i_arg+1 < argc and argv[i_arg+1][0] != '-') {
        args.push_back(argv[++i_arg]);
      }
      if (args.size() == 0) {
        std::cerr << arg << " requires arguments after it" << std::endl;
        return 1;
      } else if (args.size() == 1) {
        max_energy = std::stod(args[0]);
      } else if (args.size() == 2) {
        min_energy = std::stod(args[0]);
        max_energy = std::stod(args[1]);
      } else if (args.size() == 3) {
        min_energy = std::stod(args[0]);
        max_energy = std::stod(args[1]);
        energy_step = std::stod(args[2]);
      }
    } else if (arg == "--target") {
      std::vector<std::string> args;
      while (i_arg+1 < argc and argv[i_arg+1][0] != '-') {
        args.push_back(argv[++i_arg]);
      }
      if (args.size() != 2) {
        std::cerr << arg << " requires two arguments: Z A" << std::endl;
        return 1;
      }
      target_Z       = std::stod(args[0]);
      target_A       = std::stod(args[1]);
    } else {
      std::cout << arg << " is an unrecognized option" << std::endl;
      return 1;
    }
  }

  std::ofstream table_file(output_filename);
  if (!table_file.is_open()) {
    std::cerr << "File '" << output_filename << "' was not able to be opened." << std::endl;
    return 2;
  }

  G4double current_energy = min_energy * GeV;
  energy_step *= GeV;
  max_energy *= GeV;

  std::cout 
    << "Parameter         : Value\n"
    << "Mass A' [MeV]     : " << ap_mass*GeV << "\n"
    << "Min Energy [MeV]  : " << current_energy << "\n"
    << "Max Energy [MeV]  : " << max_energy     << "\n"
    << "Energy Step [MeV] : " << energy_step    << "\n"
    << "Lepton            : " << (muons ? "Muons" : "Electrons") << "\n"
    << "Target A [amu]    : " << target_A << "\n"
    << "Target Z [amu]    : " << target_Z << "\n"
    << std::flush;

  // the process accesses the A' mass from the G4 particle
  G4APrime::Initialize(ap_mass*GeV);
  auto model = std::make_shared<g4db::G4DarkBreMModel>("forward_only",
        0.0, 1.0, "NOT NEEDED", muons, 622, false);
  // wrap the created model in the cache so we can use it
  // to hold the xsec table and write out the CSV later
  g4db::ElementXsecCache cache(model);

  int bar_width = 80;
  int pos = 0;
  bool is_redirected = (isatty(STDOUT_FILENO) == 0);
  while (current_energy < max_energy + energy_step) {
    cache.get(current_energy, target_A, target_Z);
    current_energy += energy_step;
    if (not is_redirected) {
      int old_pos{pos};
      pos = bar_width * current_energy / max_energy;
      if (pos != old_pos) {
        std::cout << "[";
        for (int i{0}; i < bar_width; ++i) {
          if (i < pos) std::cout << "=";
          else if (i == pos) std::cout << ">";
          else std::cout << " ";
        }
        std::cout << "] " << int(current_energy / max_energy * 100.0) << " %\r";
        std::cout.flush();
      }
    }
  }
  if (not is_redirected) std::cout << std::endl;

  table_file << cache;

  table_file.close();

  return 0;
} catch (const std::exception& e) {
  std::cerr << "ERROR: " << e.what() << std::endl;
  return 127;
}
