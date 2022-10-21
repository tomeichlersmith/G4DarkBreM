/**
 * @file G4DarkBremsstrahlung.cxx
 * @brief Class providing the Dark Bremsstrahlung process class.
 * @author Michael Revering, University of Minnesota
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "G4DarkBreM/G4DarkBremsstrahlung.h"

#include "Framework/RunHeader.h"
#include "G4Electron.hh"      //for electron definition
#include "G4MuonMinus.hh"     //for muon definition
#include "G4MuonPlus.hh"      //for muon definition
#include "G4EventManager.hh"  //for EventID number
#include "G4ProcessTable.hh"  //for deactivating dark brem process
#include "G4ProcessType.hh"   //for type of process
#include "G4RunManager.hh"    //for VerboseLevel
#include "G4DarkBreM/G4DarkBreMModel.h"
#include "G4DarkBreM/DMG4Model.h"
#include "G4DarkBreM/G4APrime.h"

namespace simcore {
namespace darkbrem {

G4double ElementXsecCache::get(G4double energy, G4double A, G4double Z) {
  key_t key = computeKey(energy, A, Z);
  if (the_cache_.find(key) == the_cache_.end()) {
    if (model_.get() == nullptr) {
      EXCEPTION_RAISE("BadCache",
                      "ElementXsecCache not given a model to calculate cross "
                      "sections with.");
    }
    the_cache_[key] = model_->ComputeCrossSectionPerAtom(energy, A, Z);
  }
  return the_cache_.at(key);
}

void ElementXsecCache::stream(std::ostream& o) const {
  o << "A [au],Z [protons],Energy [MeV],Xsec [pb]\n"
    << std::setprecision(std::numeric_limits<double>::digits10 +
                         1);  // maximum precision
  for (auto const& [key, xsec] : the_cache_) {
    key_t E = key % MAX_E;
    key_t A = ((key - E) / MAX_E) % MAX_A;
    key_t Z = ((key - E) / MAX_E - A) / MAX_A;
    o << A << "," << Z << "," << E << "," << xsec / picobarn << "\n";
  }
  o << std::endl;
}

ElementXsecCache::key_t ElementXsecCache::computeKey(G4double energy,
                                                     G4double A,
                                                     G4double Z) const {
  key_t energyKey = energy;
  key_t AKey = A;
  key_t ZKey = Z;
  return (ZKey * MAX_A + AKey) * MAX_E + energyKey;
}

const std::string G4DarkBremsstrahlung::PROCESS_NAME = "eDarkBrem";
G4DarkBremsstrahlung* G4DarkBremsstrahlung::the_process_ = nullptr;

G4DarkBremsstrahlung::G4DarkBremsstrahlung(
    const framework::config::Parameters& params)
    : G4VDiscreteProcess(G4DarkBremsstrahlung::PROCESS_NAME,
                         fElectromagnetic) {
  // we need to pretend to be an EM process so the biasing framework recognizes
  // us
  the_process_ = this;
  SetProcessSubType(63);  // needs to be different from the other Em Subtypes

  global_bias_ = params.getParameter<double>("global_bias");
  only_one_per_event_ = params.getParameter<bool>("only_one_per_event");
  cache_xsec_ = params.getParameter<bool>("cache_xsec");
  ap_mass_ = params.getParameter<double>("ap_mass");
  muons_ = params.getParameter<bool>("muons");

  auto model{params.getParameter<framework::config::Parameters>("model")};
  auto model_name{model.getParameter<std::string>("name")};
  if (model_name == "vertex_library") {
    model_ = std::make_shared<G4DarkBreMModel>(model, muons_);
  } else if (model_name == "dmg4") {
    model_ = std::make_shared<DMG4Model>(model, muons_);
  } else {
    EXCEPTION_RAISE("DarkBremModel",
                    "Model named '" + model_name + "' is not known.");
  }

  // now that the model is set, calculate common xsec and put them into the
  // cache
  if (cache_xsec_) {
    element_xsec_cache_ =
        ElementXsecCache(model_);  // remake cache with model attached
    if (params.getParameter<bool>("calc_common", true)) CalculateCommonXsec();
  }
}

G4bool G4DarkBremsstrahlung::IsApplicable(const G4ParticleDefinition& p) {
  if (muons_) return &p == G4MuonMinus::Definition() or &p == G4MuonPlus::Definition();
  else return &p == G4Electron::Definition();
}

void G4DarkBremsstrahlung::PrintInfo() {
  G4cout << " Only One Per Event               : " << only_one_per_event_
         << G4endl;
  G4cout << " A' Mass [MeV]                    : " << ap_mass_ << G4endl;
  model_->PrintInfo();
}

void G4DarkBremsstrahlung::RecordConfig(ldmx::RunHeader& h) const {
  h.setIntParameter("Only One DB Per Event", only_one_per_event_);
  h.setFloatParameter("A' Mass [MeV]", ap_mass_);
  model_->RecordConfig(h);
}

G4VParticleChange* G4DarkBremsstrahlung::PostStepDoIt(const G4Track& track,
                                                       const G4Step& step) {
  // Debugging Purposes: Check if track we get is an electron
  if (not IsApplicable(*track.GetParticleDefinition()))
    EXCEPTION_RAISE(
        "DBBadTrack",
        "Dark brem process receieved a track that isn't applicable.");

  /*
   * Geant4 has decided that it is our time to interact,
   * so we are going to change the particle
   */
  ldmx_log(debug) << "A dark brem occurred!";

  if (only_one_per_event_) {
    // Deactivate the process after one dark brem if we restrict ourselves to
    // only one per event. If this is in the stepping action instead, more than
    // one brem can occur within each step. Reactivated in
    // RunManager::TerminateOneEvent Both biased and unbiased process could be
    // in the run (but not at the same time),
    //  so we turn off both while silencing the warnings from the process table.
    std::vector<G4String> db_process_name_options = {
        "biasWrapper(" + PROCESS_NAME + ")", PROCESS_NAME};
    ldmx_log(debug) << "Deactivating dark brem process";
    G4ProcessManager* pman = track.GetDefinition()->GetProcessManager();
    for (std::size_t i_proc{0}; i_proc < pman->GetProcessList()->size(); i_proc++) {
      G4VProcess* p{(*(pman->GetProcessList()))[i_proc]};
      if (p->GetProcessName().contains(PROCESS_NAME)) {
        pman->SetProcessActivation(p, false);
        break;
      }
    }
  }

  ldmx_log(debug) << "Initializing track";
  aParticleChange.Initialize(track);

  ldmx_log(debug) << "Calling model's generate change";
  model_->GenerateChange(aParticleChange, track, step);

  /*
   * Parent class has some internal counters that need to be reset,
   * so we call it before returning. It will return our shared
   * protected member variable aParticleChange that we have been modifying
   */
  ldmx_log(debug) << "Calling parent's poststepdoit";
  return G4VDiscreteProcess::PostStepDoIt(track, step);
}

void G4DarkBremsstrahlung::CalculateCommonXsec() {
  // first in pair is A, second is Z
  std::vector<std::pair<G4double, G4double>> elements = {
      std::make_pair(183.84, 74),  // tungsten
      std::make_pair(28.085, 14)   // silicon
  };

  G4double current_energy = 1.0 * MeV;
  G4double maximum_energy = 100.0 * GeV;
  if (muons_) maximum_energy *= 10;
  G4double energy_step = 1.0 * MeV;
  while (current_energy <= maximum_energy) {
    for (auto const& [A, Z] : elements)
      element_xsec_cache_.get(current_energy, A, Z);
    current_energy += energy_step;
  }
}

G4double G4DarkBremsstrahlung::GetMeanFreePath(const G4Track& track, G4double,
                                                G4ForceCondition*) {
  // won't happen if it isn't applicable
  if (not IsApplicable(*track.GetParticleDefinition())) return DBL_MAX;

  G4double energy = track.GetDynamicParticle()->GetKineticEnergy();
  G4double SIGMA = 0;
  G4Material* materialWeAreIn = track.GetMaterial();
  const G4ElementVector* theElementVector = materialWeAreIn->GetElementVector();
  const G4double* NbOfAtomsPerVolume = materialWeAreIn->GetVecNbOfAtomsPerVolume();
  
  for (size_t i = 0; i < materialWeAreIn->GetNumberOfElements(); i++) {
    G4double AtomicZ = (*theElementVector)[i]->GetZ();
    G4double AtomicA = (*theElementVector)[i]->GetA() / (g / mole);
  
    G4double element_xsec;
  
    if (cache_xsec_)
      element_xsec = element_xsec_cache_.get(energy, AtomicA, AtomicZ);
    else
      element_xsec =
          model_->ComputeCrossSectionPerAtom(energy, AtomicA, AtomicZ);
  
    SIGMA += NbOfAtomsPerVolume[i] * element_xsec;
  }
  SIGMA *= global_bias_;
  /*
  std::cout << "G4DBrem : sigma = " << SIGMA 
    << " initIntLenLeft = " << theInitialNumberOfInteractionLength
    << " nIntLenLeft = " << theNumberOfInteractionLengthLeft << std::endl;
    */
  return SIGMA > DBL_MIN ? 1. / SIGMA : DBL_MAX;
}
}  // namespace darkbrem
}  // namespace simcore
