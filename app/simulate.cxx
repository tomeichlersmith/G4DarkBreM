/**
 * @file simulate.cxx
 * definition of g4db-simulate executable
 */

#include <fstream>
#include <iostream>
#include <memory>
#include <array>

#include "QBBC.hh"
#include "G4PhysListFactory.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4UserEventAction.hh"
#include "G4UserTrackingAction.hh"
#include "G4RunManager.hh"
#include "G4MuonMinus.hh"
#include "G4Electron.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"

#include "G4DarkBreM/G4DarkBremsstrahlung.h"
#include "G4DarkBreM/G4DarkBreMModel.h"
#include "G4DarkBreM/G4APrime.h"

/**
 * example simulation application required classes
 *
 * this namespace is here _only_ to put these example classes
 * deeper into the documentation so that users can see the 
 * important stuff first.
 */
namespace g4db::example {

/**
 * basic physics constructor which simply creates the A' and the dark brem
 */
class APrimePhysics : public G4VPhysicsConstructor {
  /// handle to the process, cleaned up when the physics list is desctructed
  std::unique_ptr<G4DarkBremsstrahlung> the_process_;
  /// path to library for the model to load
  std::string library_path_;
  /// mass of A' in GeV
  double ap_mass_;
  /// true for using muons, electrons otherwise
  bool muons_;
  /// bias factor to apply everywhere
  double bias_;
 public:
  /// create the physics and store the parameters
  APrimePhysics(const std::string& lp, double m, bool mu, double b)
    : G4VPhysicsConstructor("APrime"), library_path_{lp}, ap_mass_{m}, muons_{mu}, bias_{b} {}

  /**
   * Insert A-prime into the Geant4 particle table.
   * For now we flag it as stable.
   *
   * We also define its mass here by passing the A' mass parameter.
   * Future calls to G4APrime::APrime can provide no arguments.
   *
   * Geant4 registers all instances derived from G4ParticleDefinition and
   * deletes them at the end of the run.
   */
  void ConstructParticle() final override {
    G4APrime::APrime(ap_mass_);
  }

  /**
   * Construct and configure the dark brem process
   *
   * We own the process and clean it up when the physics constructor
   * is cleaned up by Geant4 after registration.
   *
   * Lots of configuration variables here are hard-coded for this
   * simple example simulation, users of G4DarkBreM are encouraged
   * to try out the different options to see what works best for 
   * their situation.
   */
  void ConstructProcess() final override {
    the_process_ = std::unique_ptr<G4DarkBremsstrahlung>(new G4DarkBremsstrahlung(
        std::shared_ptr<g4db::G4DarkBreMModel>(new g4db::G4DarkBreMModel(
          "forward_only", /* scaling method */
          0.0, /* minimum energy threshold to dark brem [GeV] */
          1.0, /* epsilon */
          library_path_, muons_)),
        false, /* only one per event */
        bias_, /* global bias */
        true /* cache xsec */));
  }
};  // APrimePhysics

/**
 * basic 'hunk' of material in air, the material and its thickness is configurable
 *
 * The transverse (x,y) dimensions are set arbitrarily to 1m just to make
 * absolutely sure that we can contain the shower that may contain a dark brem.
 */
class Hunk : public G4VUserDetectorConstruction {
  /// depth along beam direction
  double depth_;
  /// name of material to use for volume (findable by G4NistManager)
  std::string material_;
 public:
  /**
   * Create our detector constructor, storing the configuration variables
   */
  Hunk(double d, const std::string& m)
    : G4VUserDetectorConstruction(), depth_{d}, material_{m} {}

  /**
   * Construct the geometry
   *
   * We build the world only slighly larger than the single hunk
   * of material at its center. The hunk is shifted to be
   * downstream (along z) of the origin so that the primary generator
   * can simply shoot from the origin along z.
   */
  virtual G4VPhysicalVolume* Construct() final override {
    // Get nist material manager
    G4NistManager* nist = G4NistManager::Instance();
    using CLHEP::mm;

    G4double box_half_x{500*mm},
             box_half_y{500*mm},
             box_half_z{depth_/2*mm};
    G4Material* box_mat = nist->FindOrBuildMaterial(material_);
    if (not box_mat) {
      throw std::runtime_error("Material '"+material_+"' unknown to G4NistManager.");
    }

    G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
    if (not world_mat) {
      throw std::runtime_error("Material 'G4_AIR' unknown to G4NistManager.");
    }

    G4double world_half_z = 2*box_half_z+2;
    G4Box* solidWorld =
      new G4Box("World", 1.1*box_half_x, 1.1*box_half_y, world_half_z);

    G4LogicalVolume* logicWorld =
      new G4LogicalVolume(solidWorld,
          world_mat,
          "World");

    G4VPhysicalVolume* physWorld =
      new G4PVPlacement(0, //no rotation
          G4ThreeVector(), // center nudged upstream a bit
          logicWorld,      //its logical volume
          "World",         //its name
          0,               //its mother  volume
          false,           //no boolean operation
          0,               //copy number
          false);          //overlaps checking

    G4Box* solidBox = new G4Box("Box",
        box_half_x, box_half_y, box_half_z);

    G4LogicalVolume* logicBox = new G4LogicalVolume(solidBox,
        box_mat, "Box");

    // providing mother volume attaches us to the world volume
    new G4PVPlacement(0, //no rotation
        G4ThreeVector(0.,0.,box_half_z+1), //at (0,0,box_half_z+1)
        logicBox,        //its logical volume
        "Envelope",      //its name
        logicWorld,      //its mother  volume
        false,           //no boolean operation
        0,               //copy number
        false);          //overlaps checking

    //always return the physical World
    return physWorld;
  }
};

/**
 * The event information we care about for studying the model
 *
 * Since a new instance of this object is created for each
 * event (and destroyed at the end of the event). The default
 * values of its members correspond to the starting values
 * at the beginning of an event.
 */
class OutgoingKinematics : public G4VUserEventInformation {
  /// have we found the dark brem products yet?
  bool found_{false};
  /// the four momentum of the recoil lepton
  std::array<double,4> recoil_;
  /// the four momentum of the produced dark photon (A')
  std::array<double,4> aprime_;
 public:
  /**
   * Check if the dark brem has been found
   */
  bool found() const {
    return found_;
  }

  /**
   * Set the recoil 4 momentum based on the passed track
   */
  void setRecoil(const G4Track* track) {
    found_ = true;
    recoil_ = {
      track->GetTotalEnergy(),
      track->GetMomentum().x(),
      track->GetMomentum().y(),
      track->GetMomentum().z()
    };
  }

  /**
   * Set the dark photon 4 momentum based on the passed track
   */
  void setAPrime(const G4Track* track) {
    found_ = true;
    aprime_ = {
      track->GetTotalEnergy(),
      track->GetMomentum().x(),
      track->GetMomentum().y(),
      track->GetMomentum().z()
    };
  }

  /**
   * Write out the two four-momenta in CSV format to the input stream
   */
  void stream(std::ostream& o) const {
    o << recoil_[0] << ',' << recoil_[1] << ',' << recoil_[2] << ',' << recoil_[3] << ','
      << aprime_[0] << ',' << aprime_[1] << ',' << aprime_[2] << ',' << aprime_[3] << '\n';
  }
  
  /**
   * Reequired by Geant4, simply print to std::cout using stream
   */
  virtual void Print() const final override {
    stream(std::cout);
    std::cout << std::flush;
  }
  
  /**
   * Helper function to retrieve the current instance of the event information
   */
  static OutgoingKinematics* get() {
    auto user_info = G4EventManager::GetEventManager()->GetUserInformation();
    if (not user_info) throw std::runtime_error("Attempting to 'get' not-created event information");
    return dynamic_cast<OutgoingKinematics*>(user_info);
  }
};  // OutgoingKinematics

/**
 * the primary generator, a simple particle gun restricted to electrons or muons
 * along the z axis
 */
class LeptonBeam : public G4VUserPrimaryGeneratorAction {
  /// the gun we use for the beam
  G4ParticleGun gun_;
 public:
  /**
   * Configure the beam to be of the input energy and lepton
   *
   * Shoot along the z axis, the energy is in GeV and we
   * shoot from the origin.
   */
  LeptonBeam(double e, bool muons)
    : G4VUserPrimaryGeneratorAction() {
      if (muons) gun_.SetParticleDefinition(G4MuonMinus::MuonMinus());
      else gun_.SetParticleDefinition(G4Electron::Electron());
      gun_.SetParticleEnergy(e*GeV);
      gun_.SetParticlePosition(G4ThreeVector());
      gun_.SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
    }

  /**
   * Start an event by providing primaries
   *
   * We also construct the OutgoingKinematics structure for this event
   */
  void GeneratePrimaries(G4Event* event) final override {
    event->SetUserInformation(new OutgoingKinematics);
    gun_.GeneratePrimaryVertex(event);
  }
};

/**
 * event action used to store the OutgoingKinematics *if* a dark brem occurred
 *
 * This uses OutgoingKinematics::stream to write out the CSV where data will be stored.
 * 
 * @note The stream method at the CSV header row written here need to match for the CSV
 * to make sense.
 *
 * We don't do any caching, just trusting the std::ofstream to handle the caching,
 * only flushing when necessary and when destructed.
 *
 * We also print out the number of events that successfully had a dark brem compared
 * to the number of events requested. This is helpful for the user so that they
 * know (1) there is not a problem and (2) potential tuning of the bias factor.
 */
class PersistDarkBremProducts : public G4UserEventAction {
  /// the output file we are writing to
  std::ofstream out_;
  /// number of events that we simulated
  long unsigned int events_started_{0};
  /// number of events with a dark brem in it
  long unsigned int events_completed_{0};
 public:
  /**
   * Open the output CSV and write the header row
   */
  PersistDarkBremProducts(const std::string& out_file)
    : G4UserEventAction(), out_{out_file} {
    if (not out_.is_open()) {
      throw std::runtime_error("Unable to open output file '"+out_file+"'.");
    }
    out_ << "recoil_energy,recoil_px,recoil_py,recoil_pz,aprime_energy,aprime_px,aprime_py,aprime_pz\n";
    out_.flush();
  }

  /**
   * Print out the number of events with a dark brem compared to the requested number
   */
  ~PersistDarkBremProducts() {
    std::cout << "[g4db-simulate] Able to generate a dark brem " 
      << events_completed_ << " / " << events_started_ 
      << " events" << std::endl;
  }
  
  /**
   * Check the OutgoingKinematics and write out the four-momenta if the dark brem ocurred
   */
  void EndOfEventAction(const G4Event* event) final override {
    ++events_started_;
    auto ek{dynamic_cast<OutgoingKinematics*>(event->GetUserInformation())};
    if (ek->found()) {
      ++events_completed_;
      ek->stream(out_);
    }
  }
};  // PersistDarkBremProducts

/**
 * Look through the tracks to find the dark brem products
 */
class FindDarkBremProducts : public G4UserTrackingAction {
 public:
  /**
   * Check the new track that has been provided to see if it was
   * created by the dark brem process. If it was, pass this track
   * onto OutgoingKinematics depending on if it is the A' or the lepton.
   */
  void PreUserTrackingAction(const G4Track* track) final override {
    const G4VProcess* creator = track->GetCreatorProcess();
    if (creator and creator->GetProcessName().contains(
                        G4DarkBremsstrahlung::PROCESS_NAME)) {
      if (track->GetParticleDefinition() == G4APrime::APrime()) {
        OutgoingKinematics::get()->setAPrime(track);
      } else {
        OutgoingKinematics::get()->setRecoil(track);
      }
    } // track created by dark brem process
  }
};  // FindDarkBremProducts

}  // namespace g4db::example

/**
 * print out how to use g4db-simulate
 */
void usage() {
  std::cout <<
    "\n"
    "USAGE\n"
    "  g4db-simulate [options] DB-LIB NUM-EVENTS\n"
    "\n"
    "ARGUMENTS\n"
    "  DB-LIB     : dark brem library to scale from\n"
    "               the user is expected to make sure that this argument aligns with\n"
    "               the other options (lepton, incident beam energy, A' mass, etc...)\n"
    "  NUM-EVENTS : number of events to **request**\n"
    "               since Geant4 decides when a dark brem will occurr, it is important\n"
    "               to allow some beam leptons to /not/ dark brem in the target so a realistic\n"
    "               distribution of dark brem vertices is sampled."
    "\n"
    "OPTIONS\n"
    "  -h, --help    : print this usage and exit"
    "  --muons       : run using muons (without this flag, assumes electrons)\n"
    "  -m, --ap-mass : mass of the dark photon (A') in GeV (defaults to 0.1 for electrons and 1. for muons)\n"
    "  -d, --depth   : thickness of target in mm (defaults to 18 for electrons, 2000 for muons)\n"
    "  -t, --target  : target material, must be findable by G4NistManager\n"
    "                  (defaults to G4_W for electrons and G4_Cu for muons)\n"
    "  -o, --output  : output file to write CSV data to (defaults to 'events.csv')\n"
    "  -b, --bias    : biasing factor to use to encourage dark brem\n"
    "                  a good starting point is generally the A' mass squared, so that is the default\n"
    "  -e, --beam    : Beam energy in GeV (defaults to 4 for electrons and 100 for muons)\n"
    "  --mat-list    : print the full list from G4NistManager and exit\n"
    "\n"
    << std::flush;
}

/**
 * definition of g4db-simulate
 *
 * After parsing the command line arguments, we simply do the 
 * standard initialization and running procedure for Geant4.
 */
int main(int argc, char* argv[]) try {
  bool muons{false};
  double depth{-1};
  std::string target{};
  std::string output{"events.csv"};
  double bias{-1.};
  double beam{-1.};
  double ap_mass{-1.};
  std::vector<std::string> positional;
  for (int i_arg{1}; i_arg < argc; ++i_arg) {
    std::string arg{argv[i_arg]};
    if (arg == "-h" or arg == "--help") {
      usage();
      return 0;
    } else if (arg == "--mat-list") {
      auto nist{G4NistManager::Instance()};
      nist->ListMaterials("simple");
      nist->ListMaterials("compound");
      nist->ListMaterials("hep");
      return 0;
    } else if (arg == "--muons") {
      muons = true;
    } else if (arg == "-o" or arg == "--output") {
      if (i_arg+1 >= argc) {
        std::cerr << arg << " requires an argument after it" << std::endl;
        return 1;
      }
      output = argv[++i_arg];
    } else if (arg == "-t" or arg == "--target") {
      if (i_arg+1 >= argc) {
        std::cerr << arg << " requires an argument after it" << std::endl;
        return 1;
      }
      target = argv[++i_arg];
    } else if (arg == "-m" or arg == "--ap-mass") {
      if (i_arg+1 >= argc) {
        std::cerr << arg << " requires an argument after it" << std::endl;
        return 1;
      }
      ap_mass = std::stod(argv[++i_arg]);
    } else if (arg == "-d" or arg == "--depth") {
      if (i_arg+1 >= argc) {
        std::cerr << arg << " requires an argument after it" << std::endl;
        return 1;
      }
      depth = std::stod(argv[++i_arg]);
    } else if (arg == "-b" or arg == "--bias") {
      if (i_arg+1 >= argc) {
        std::cerr << arg << " requires an argument after it" << std::endl;
        return 1;
      }
      bias = std::stod(argv[++i_arg]);
    } else if (arg == "-e" or arg == "--beam") {
      if (i_arg+1 >= argc) {
        std::cerr << arg << " requires an argument after it" << std::endl;
        return 1;
      }
      beam = std::stod(argv[++i_arg]);
    } else if (arg[0] == '-') {
      std::cerr << arg << " is not a recognized option" << std::endl;
      return 1;
    } else {
      positional.push_back(arg);
    }
  }

  if (positional.size() != 2) {
    std::cerr << "Exactly two positional arguments are required: DB-LIB NUM-EVENTS" << std::endl;
    return 1;
  }

  std::string db_lib = positional[0];
  int num_events = std::stoi(positional[1]);

  if (muons) {
    if (ap_mass < 0.) ap_mass = 1.;
    if (beam < 0.) beam = 100.;
    if (depth < 0.) depth = 2000;
    if (target.empty()) target = "G4_Cu";
  } else {
    if (ap_mass < 0.) ap_mass = 0.1;
    if (beam < 0.) beam = 4.;
    if (depth < 0.) depth = 18;
    if (target.empty()) target = "G4_W";
  }

  if (bias < 0.) bias = ap_mass*ap_mass;

  auto run = std::unique_ptr<G4RunManager>(new G4RunManager);

  run->SetUserInitialization(new g4db::example::Hunk(depth,target));

  G4VModularPhysicsList* physics = new QBBC;
  physics->RegisterPhysics(new g4db::example::APrimePhysics(db_lib, ap_mass, muons, bias));
  run->SetUserInitialization(physics);

  run->Initialize();

  run->SetUserAction(new g4db::example::FindDarkBremProducts);
  run->SetUserAction(new g4db::example::PersistDarkBremProducts(output));
  run->SetUserAction(new g4db::example::LeptonBeam(beam, muons));

  run->BeamOn(num_events);

  return 0;
} catch (const std::exception& e) {
  std::cerr << "ERROR: " << e.what() << std::endl;
  return 127;
}
