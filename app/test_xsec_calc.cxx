
//----------------//
//   C++ StdLib   //
//----------------//
#include <fstream>
#include <iostream>

#include <boost/program_options.hpp>

//----------//
//   LDMX   //
//----------//
#include "Framework/Configure/Parameters.h"
#include "G4DarkBreM/G4DarkBremsstrahlung.h"
#include "G4DarkBreM/G4APrime.h"

template<typename T>
const T& get(boost::program_options::variables_map& vm, const std::string& name, const T& def) {
  if (vm.count(name)) return vm[name].as<T>();
  else return def;
}

/**
 * The executable main for printing out the table.
 */
int main(int argc, char* argv[]) try {
  boost::program_options::options_description desc(
      "Calculate dark brem cross sections and write them out to a CSV table"
      );
  desc.add_options()
    ("help,h", 
      "produce this help and exit")
    ("output,o", 
      boost::program_options::value<std::string>()->default_value("xsec.csv"), 
      "output file to write xsec to")
    ("model,m", 
      boost::program_options::value<std::string>()->default_value("g4db"), 
      "model to use for calculating")
    ("ap-mass",
      boost::program_options::value<double>(),
      "mass of A' in MeV (defaults to 100 for electrons and 200 for muons)")
    ("muons",
      "use muons as incident lepton rather than electrons")
    ("min-energy",
      boost::program_options::value<double>(),
      "minimum energy [GeV] to start xsec calculation from (defaults to 2*ap-mass)")
    ("max-energy",
      boost::program_options::value<double>(),
      "maximum energy [GeV] to calculate xsec to (defaults to 100GeV for electrons and 1000GeV for muons)")
    ("energy-step",
      boost::program_options::value<double>()->default_value(1000.),
      "difference between adjacent steps in energy [MeV]")
    ("target-density",
      boost::program_options::value<double>(),
      "Density of target nucleus (defaults to 19.3 for electrons and 8.96 for muons)")
    ("target-Z",
      boost::program_options::value<int>(),
      "Z of target nucleus (defaults to 74 for electrons and 29 for muons)")
    ("target-A",
      boost::program_options::value<double>(),
      "A of target nucleus (defaults to 183.84 for electrons and 64.0 for muons)")
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

  std::string model_name{vm["model"].as<std::string>()};
  if (model_name == "g4db") {
    model_name = "vertex_library";
  } else if (model_name != "dmg4") {
    std::cerr << "Model '" << model_name << "' is not 'g4db' or 'dmg4'." << std::endl;
    return -1;
  }

  std::ofstream table_file(vm["output"].as<std::string>());
  if (!table_file.is_open()) {
    std::cerr << "File '" << vm["output"].as<std::string>() << "' was not able to be opened."
              << std::endl;
    return 2;
  }

  bool muons = vm.count("muons") > 0;

  G4double ap_mass, max_energy;
  double target_A, density;
  int target_Z;
  if (muons) {
    ap_mass     = get(vm, "ap-mass"   , 200. ) * MeV;
    max_energy  = get(vm, "max-energy", 1000.) * GeV;
    density     = get(vm, "density"   , 8.96 );
    target_Z    = get(vm, "target-Z"  , 29   );
    target_A    = get(vm, "target-A"  , 64.0 );
  } else {
    ap_mass     = get(vm, "ap-mass"   , 100.   ) * MeV;
    max_energy  = get(vm, "max-energy", 100.   ) * GeV;
    density     = get(vm, "density"   , 19.30  );
    target_Z    = get(vm, "target-Z"  , 74     );
    target_A    = get(vm, "target-A"  , 183.84 );
  }

  framework::config::Parameters model;
  model.addParameter<std::string>("name", model_name);
  model.addParameter("epsilon", 1.);
  if (model_name == "vertex_library") {
    model.addParameter<std::string>("library_path", "NOTNEEDED");
    model.addParameter<std::string>("method", "forward_only");
    model.addParameter("threshold", 0.0);
    model.addParameter("load_library", false);
  } else {
    // DMG4 requires target-material at construction time
    model.addParameter("density", density);
    model.addParameter("targetA", target_A);
    model.addParameter<double>("targetZ", target_Z);
  }

  framework::config::Parameters process;
  process.addParameter("model", model);
  process.addParameter("ap_mass", ap_mass/MeV);
  process.addParameter("enable", true);
  process.addParameter("only_one_per_event", false);
  process.addParameter("cache_xsec", true);
  process.addParameter("global_bias", 1.);
  process.addParameter("muons", muons);
  process.addParameter("calc_common",false);

  G4double current_energy = get(vm, "min-energy", 2 * ap_mass / GeV) * GeV;
  G4double energy_step = vm["energy-step"].as<double>() * MeV;

  // the process accesses the A' mass from the G4 particle
  simcore::darkbrem::G4APrime::APrime(process.getParameter<double>("ap_mass"));
  // create the process to do proper initializations
  //    this calculates "common" cross sections as well
  simcore::darkbrem::G4DarkBremsstrahlung db_process(process);

  int bar_width = 80;
  int pos = 0;
  bool is_redirected = (isatty(STDOUT_FILENO) == 0);
  while (current_energy < max_energy + energy_step) {
    db_process.getCache().get(current_energy, target_A, target_Z);
    current_energy += energy_step;
    int old_pos{pos};
    pos = bar_width * current_energy / max_energy;
    if (pos != old_pos) {
      std::cout << "[";
      for (int i{0}; i < bar_width; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
      }
      std::cout << "] " << int(current_energy / max_energy * 100.0) << " %";
      if (is_redirected) std::cout << "\n";
      else std::cout << "\r";
      std::cout.flush();
    }
  }
  std::cout << std::endl;

  table_file << db_process.getCache();

  table_file.close();

  db_process.PrintInfo();

  return 0;
} catch (const std::exception& e) {
  std::cerr << "ERROR: " << e.what() << std::endl;
  return 127;
}
