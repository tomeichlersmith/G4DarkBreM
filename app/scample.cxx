
#include <iostream>
#include <fstream>

#include <boost/program_options.hpp>


#include "G4Electron.hh"
#include "G4MuonMinus.hh"

#include "G4DarkBreM/G4DarkBreMModel.h"
#include "G4DarkBreM/G4APrime.h"

template<typename T>
const T& get(boost::program_options::variables_map& vm, const std::string& name, const T& def) {
  if (vm.count(name)) return vm[name].as<T>();
  else return def;
}

int main(int argc, char* argv[]) try {
  boost::program_options::options_description desc(
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
      "OPTIONS"
      );
  desc.add_options()
    ("help,h", "produce this help and exit")
    ("output,o", boost::program_options::value<std::string>()->default_value("scaled.root"),
     "output file to write scaled outgoing kinematics to")
    ("incident-energy,E", boost::program_options::value<double>(),
     "incident energy in GeV (defaults to 4 for electrons and 100 for muons)")
    ("num-events,N", boost::program_options::value<int>()->default_value(100),
     "number of dark brems to simulate")
    ("db-lib,L", boost::program_options::value<std::string>(),
     "DB event library to load")
    ("ap-mass",
      boost::program_options::value<double>(),
      "mass of A' in MeV (defaults to 100 for electrons and 1000 for muons)")
    ("muons",
      "use muons as incident lepton rather than electrons")
  ;

  boost::program_options::variables_map vm;
  boost::program_options::store(
      boost::program_options::command_line_parser(argc, argv)
      .options(desc).run(), vm);
  boost::program_options::notify(vm);

  if (vm.count("help")) {
    std::cout << desc;
    return 0;
  }

  if (vm.count("db-lib") == 0) {
    std::cerr << "ERROR: DB event library not provided." << std::endl;
    return 1;
  }

  int n_events = vm["num-events"].as<int>();
  bool muons = vm.count("muons") > 0;
  double ap_mass, incident_energy, lepton_mass;
  if (muons) {
    ap_mass     = get(vm, "ap-mass"    , 1000. ) * MeV;
    incident_energy = get(vm, "incident-energy", 100. );
    lepton_mass = G4MuonMinus::MuonMinus()->GetPDGMass() / GeV;
  } else {
    ap_mass     = get(vm, "ap-mass"    , 100. ) * MeV;
    incident_energy = get(vm, "incident-energy", 4.   );
    lepton_mass = G4Electron::Electron()->GetPDGMass() / GeV;
  }

  // the process accesses the A' mass from the G4 particle
  G4APrime::APrime(ap_mass/MeV);
  // create the model, this is where the LHE file is parsed
  //    into an in-memory library to sample and scale from
  g4db::G4DarkBreMModel db_model("forward_only", 0.0, 1.0, 
      vm["db-lib"].as<std::string>(), muons);
  db_model.PrintInfo();
  printf("   %-16s %f\n", "Lepton Mass [MeV]:", lepton_mass);
  printf("   %-16s %f\n", "A' Mass [MeV]:", ap_mass/MeV);

  std::ofstream f{vm["output"].as<std::string>()};
  if (not f.is_open()) {
    std::cerr << "Unable to open output file for writing." << std::endl;
    return -1;
  }
  f << "recoil_energy,recoil_px,recoil_py,recoil_pz\n";

  for (int i_event{0}; i_event < n_events; ++i_event) {
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
