/** 
 * @file scample.cxx
 * definition of g4db-scample executable
 */

#include <iostream>
#include <fstream>

#include "G4Electron.hh"
#include "G4MuonMinus.hh"

#include "G4DarkBreM/G4DarkBreMModel.h"
#include "G4DarkBreM/G4APrime.h"

/**
 * printout how to use g4db-scample
 */
void usage() {
  std::cout << 
      "USAGE:\n"
      "  g4db-scample [options] db-lib\n"
      "\n"
      "Run the scaling procedure for the input incident energy and madgraph file\n"
      "\n"
      "This executable is a low-level way to directly test the scaling procedure implemented\n"
      "inside the G4DarkBreMModel without cluttering the results with the rest of the Geant4\n"
      "simulation machinery. This means a better understanding of how the model functions is\n"
      "necessary to be able to effectively use this program.\n"
      " - The 'incident energy' input here is the energy of the lepton JUST BEFORE it dark brems.\n"
      " - The scaling procedure should scale from a MG sample at an energy ABOVE the incident energy\n"
      " - The scaling procedure generates the recoil lepton's kinematics assuming the incident\n"
      "   lepton is traveling along the z-axis. The user is expected to rotate to the actual incident\n"
      "   frame and calculate the outgoing dark photon kinematics assuming conservation of momentum.\n"
      "\n"
      "ARGUMENTS\n"
      "  db-lib : dark brem event library to load and sample\n"
      "\n"
      "OPTIONS\n"
      "  -h,--help             : produce this help and exit\n"
      "  -o,--output           : output file to write scaled events to\n"
      "  -E,--incident-energy  : energy of incident lepton in GeV\n"
      "  -N,--num-events       : number of events to sample and scale\n"
      "  -M,--ap-mass          : mass of dark photon in GeV\n"
      "  --muons               : pass to set lepton to muons (otherwise electrons)\n"
      << std::flush;
}

/**
 * definition of g4db-scample
 *
 * We only need to configure the G4DarkBreMModel so
 * we simply define G4APrime and then construct the model
 * so we can call G4DarkBreMModel::scample for the input
 * number of events.
 */
int main(int argc, char* argv[]) try {
  std::string output_filename{"scaled.csv"};
  double incident_energy{4};
  int num_events{10};
  std::string db_lib;
  double ap_mass{0.1};
  bool muons{false};
  for (int i_arg{1}; i_arg < argc; i_arg++) {
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
    } else if (arg == "-E" or arg == "--incident-energy") {
      if (i_arg+1 >= argc) {
        std::cerr << arg << " requires an argument after it" << std::endl;
        return 1;
      }
      incident_energy = std::stod(argv[++i_arg]);
    } else if (arg == "-M" or arg == "--ap-mass") {
      if (i_arg+1 >= argc) {
        std::cerr << arg << " requires an argument after it" << std::endl;
        return 1;
      }
      ap_mass = std::stod(argv[++i_arg]);
    } else if (arg == "-N" or arg == "--num-events") {
      if (i_arg+1 >= argc) {
        std::cerr << arg << " requires an argument after it" << std::endl;
        return 1;
      }
      num_events = std::stoi(argv[++i_arg]);
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

  double lepton_mass;
  if (muons) {
    lepton_mass = G4MuonMinus::MuonMinus()->GetPDGMass() / GeV;
  } else {
    lepton_mass = G4Electron::Electron()->GetPDGMass() / GeV;
  }

  // the process accesses the A' mass from the G4 particle
  G4APrime::APrime(ap_mass/GeV);
  // create the model, this is where the LHE file is parsed
  //    into an in-memory library to sample and scale from
  g4db::G4DarkBreMModel db_model("forward_only", 0.0, 1.0, db_lib, muons);
  db_model.PrintInfo();
  printf("   %-16s %f\n", "Lepton Mass [MeV]:", lepton_mass);
  printf("   %-16s %f\n", "A' Mass [MeV]:", ap_mass/MeV);

  std::ofstream f{output_filename};
  if (not f.is_open()) {
    std::cerr << "Unable to open output file for writing." << std::endl;
    return -1;
  }
  f << "recoil_energy,recoil_px,recoil_py,recoil_pz\n";

  for (int i_event{0}; i_event < num_events; ++i_event) {
    G4ThreeVector recoil = db_model.scample(incident_energy, lepton_mass);
    double recoil_energy = sqrt(recoil.mag2() + lepton_mass*lepton_mass);

    f << recoil_energy << ','
      << recoil.x() << ','
      << recoil.y() << ','
      << recoil.z() << '\n';
  }

  f.close();
  return 0;
} catch (const std::exception& e) {
  std::cerr << "ERROR: " << e.what() << std::endl;
  return 127;
}
