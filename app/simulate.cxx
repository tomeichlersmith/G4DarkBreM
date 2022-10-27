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
 * basic physics constructor which simply creates the A' and the dark brem
 */
class APrimePhysics : public G4VPhysicsConstructor {
  std::unique_ptr<G4DarkBremsstrahlung> the_process_;
  std::string library_path_;
  double ap_mass_;
  bool muons_;
  double bias_;
 public:
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
        muons_,
        false, /* only one per event */
        bias_, /* global bias */
        true /* cache xsec */));
    the_process_->SetModel(
        std::shared_ptr<g4db::G4DarkBreMModel>(new g4db::G4DarkBreMModel(
          "forward_only", /* scaling method */
          0.0, /* minimum energy threshold to dark brem [GeV] */
          1.0, /* epsilon */
          library_path_, muons_)));
  }
};  // APrimePhysics

/**
 * basic 'hunk' of material in air, the material and its thickness is configurable
 *
 * The transverse (x,y) dimensions are set arbitrarily to 1m just to make
 * absolutely sure that we can contain the shower that may contain a dark brem.
 */
class Hunk : public G4VUserDetectorConstruction {
  double depth_;
  std::string material_;
 public:
  Hunk(double d, const std::string& m)
    : G4VUserDetectorConstruction(), depth_{d}, material_{m} {}
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

class OutgoingKinematics : public G4VUserEventInformation {
  bool found_{false};
  std::array<double,4> recoil_;
  std::array<double,4> aprime_;
 public:
  OutgoingKinematics()
    : G4VUserEventInformation() {}
  ~OutgoingKinematics() = default;
  bool found() const {
    return found_;
  }
  void setRecoil(const G4Track* track) {
    found_ = true;
    recoil_ = {
      track->GetTotalEnergy(),
      track->GetMomentum().x(),
      track->GetMomentum().y(),
      track->GetMomentum().z()
    };
  }
  void setAPrime(const G4Track* track) {
    found_ = true;
    aprime_ = {
      track->GetTotalEnergy(),
      track->GetMomentum().x(),
      track->GetMomentum().y(),
      track->GetMomentum().z()
    };
  }
  void stream(std::ostream& o) const {
    o << recoil_[0] << ',' << recoil_[1] << ',' << recoil_[2] << ',' << recoil_[3] << ','
      << aprime_[0] << ',' << aprime_[1] << ',' << aprime_[2] << ',' << aprime_[3] << '\n';
  }
  virtual void Print() const final override {
    stream(std::cout);
    std::cout << std::flush;
  }
  static OutgoingKinematics* get() {
    auto user_info = G4EventManager::GetEventManager()->GetUserInformation();
    if (not user_info) throw std::runtime_error("Attempting to 'get' not-created event information");
    return dynamic_cast<OutgoingKinematics*>(user_info);
  }
};

/**
 * the primary generator, a simple particle gun restricted to electrons or muons
 * along the z axis
 */
class LeptonBeam : public G4VUserPrimaryGeneratorAction {
  G4ParticleGun gun_;
 public:
  LeptonBeam(double e, bool muons)
    : G4VUserPrimaryGeneratorAction() {
      if (muons) gun_.SetParticleDefinition(G4MuonMinus::MuonMinus());
      else gun_.SetParticleDefinition(G4Electron::Electron());
      gun_.SetParticleEnergy(e*GeV);
      gun_.SetParticlePosition(G4ThreeVector());
      gun_.SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
    }
  void GeneratePrimaries(G4Event* event) final override {
    event->SetUserInformation(new OutgoingKinematics);
    gun_.GeneratePrimaryVertex(event);
  }
};

/**
 * our tracking+event action so we can save the outgoing kinematics of the dark brem
 */
class PersistDarkBremProducts : public G4UserEventAction {
  std::ofstream out_;
 public:
  PersistDarkBremProducts(const std::string& out_file)
    : G4UserEventAction(), out_{out_file} {
    if (not out_.is_open()) {
      throw std::runtime_error("Unable to open output file '"+out_file+"'.");
    }
    out_ << "recoil_energy,recoil_px,recoil_py,recoil_pz,aprime_energy,aprime_px,aprime_py,aprime_pz\n";
    out_.flush();
  }
  void EndOfEventAction(const G4Event* event) final override {
    auto ek{dynamic_cast<OutgoingKinematics*>(event->GetUserInformation())};
    if (ek->found()) ek->stream(out_);
  }
};  // PersistDarkBremProducts

class FindDarkBremProducts : public G4UserTrackingAction {
 public:
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

/**
 * The executable main for printing out the table.
 */
int main(int argc, char* argv[]) try {
  auto run = std::unique_ptr<G4RunManager>(new G4RunManager);

  run->SetUserInitialization(new Hunk(18.,"G4_W"));

  //G4PhysListFactory factory;
  G4VModularPhysicsList* physics = new QBBC; //factory.GetReferencePhysList("FTFP_BERT");
  physics->RegisterPhysics(new APrimePhysics(argv[1], 100., false, 10000.));
  run->SetUserInitialization(physics);

  run->Initialize();

  run->SetUserAction(new FindDarkBremProducts);
  run->SetUserAction(new PersistDarkBremProducts("events.csv"));
  run->SetUserAction(new LeptonBeam(4.0, false));

  run->BeamOn(10);

  return 0;
} catch (const std::exception& e) {
  std::cerr << "ERROR: " << e.what() << std::endl;
  return 127;
}
