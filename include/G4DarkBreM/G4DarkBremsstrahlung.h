/**
 * @file G4DarkBremsstrahlung.h
 * @brief Class providing the Dark Bremsstrahlung process class.
 * @author Michael Revering, University of Minnesota
 * @author Tom Eichlersmith, University of Minnesota
 */

#ifndef SIMCORE_DARKBREM_G4DARKBREMSSTRAHLUNG_H_
#define SIMCORE_DARKBREM_G4DARKBREMSSTRAHLUNG_H_

#include "G4DarkBreM/Parameters.h"

// Geant
#include "G4VDiscreteProcess.hh"

class G4String;
class G4ParticleDefinition;

namespace g4db {

/**
 * @class G4DarkBremsstrahlung
 *
 * Class that represents the dark brem process.
 * An electron is allowed to brem a dark photon
 *
 * @TODO allow positrons to dark brem as well
 */
class G4DarkBremsstrahlung : public G4VDiscreteProcess {
 public:
  /**
   * The name of this process in Geant4
   *
   * @note This process name should be used in all places that
   * can depend on this file. (For example RunManager and DetectorConstruction).
   * There are other places that _can't_ use this constant
   * directly, so if this name is changed you also need to change
   * the following places.
   *  - Python: Biasing.filters.TrackProcessFilter.dark_brem
   *  - C++: Event/SimParticle::createProcessMap
   */
  static const std::string PROCESS_NAME;

  /**
   * The created instance. This will exist if dark brem has been enabled and _after_
   * run initialization where the physics lists are constructred
   */
  static G4DarkBremsstrahlung* Get() { return the_process_; }

  /**
   * Constructor
   *
   * Configures this process by doing three main things:
   *  1. Registers this process with Geant4 as a 'fElectromagnetic' process
   *     - Needed for Geant4 Biasing Framework to recognize this process as
   * "bias-able"
   *  2. Defines the EM subtype as one different from all other EM processes
   *     - Needed so we don't replace another EM process
   *  3. Configures the process and passes the model parameters to the model
   *
   * If caching the cross section is enabled, we calculate several
   * common cross sections immediately to help even out the time
   * it takes to simulate events.
   * @see CalculateCommonXsec
   */
  G4DarkBremsstrahlung(const framework::config::Parameters& params);

  /**
   * Destructor
   */
  virtual ~G4DarkBremsstrahlung() = default;

  /**
   * Checks if the passed particle should be able to do this process
   *
   * @return true if particle is electron
   */
  virtual G4bool IsApplicable(const G4ParticleDefinition& p);

  /**
   * Reports the parameters to G4cout.
   *
   * @see G4DarkBremsstrahlungModel::PrintInfo
   */
  virtual void PrintInfo();

  /**
   * This is the function actually called by Geant4 that does the dark brem
   * interaction.
   *
   * aParticleChange is a protected member variable of G4VDiscreteProcess
   * that we should edit here.
   *
   * If only one per event is set, then we deactivate the dark brem process,
   * ensuring only one dark brem per step and per event.
   * Reactivated in RunManager::TerminateOneEvent.
   *
   * @see RunManager::TerminateOneEvent
   * @see G4DarkBremsstrahlungModel::GenerateChange
   * @param[in] track current G4Track that is being stepped
   * @param[in] step current step that just finished
   * @returns G4VParticleChange detailing how this process changes the track
   */
  virtual G4VParticleChange* PostStepDoIt(const G4Track& track,
                                          const G4Step& step);

  /**
   * Calculate common cross sections for the cache using the already-created
   * model.
   *
   * This method is public so that we can access it for writing a short
   * executable to print the cache to a file.
   */
  void CalculateCommonXsec();

  /**
   * Get a reference to the cross section cache.
   *
   * Again, this method is public only to be available to the executable
   * that generates a cross section table and testing.
   * Do not use this unless you really know what you are doing.
   */
  ElementXsecCache& getCache() { return element_xsec_cache_; }

 protected:
  /**
   * Calculate the mean free path given the input conditions
   *
   * We maintain a cache for the cross sections calculated by the model
   * so that later in the run it is less likely that the model will
   * need to be called to calculate the cross section. This is done
   * in order to attempt to improve speed of simulation and avoid
   * repetition of the same, deterministic calculations.
   *
   * If you want to turn off the cache-ing behavior, set 'cache_xsec' to false
   * in the python configuratin for the dark brem process.
   *
   * @see G4DarkBremsstrahlungModel::ComputeCrossSectionPerAtom
   * @param[in] track G4Track that is being stepped
   * @param[in] prevStepSize G4double measuring previous step size, unused
   * @param[in] condition G4ForceCondition, always NotForced for
   * G4VDiscreteProcess, unused
   * @returns G4double mean free path of the particle
   */
  G4double GetMeanFreePath(const G4Track& track, G4double prevStepSize,
                           G4ForceCondition* condition);

 private:
  /** remove ability to assign this object */
  G4DarkBremsstrahlung& operator=(const G4DarkBremsstrahlung& right);

  /** remove ability to copy construct */
  G4DarkBremsstrahlung(const G4DarkBremsstrahlung&);

  /**
   * Only allow the dark brem to happen once per event.
   *
   * This allows for the dark brem process to be de-activated when
   * SampleSecondaries is called.
   *
   * The dark brem process is _always_ re-activated in the
   * RunManager::TerminateOneEvent method. This reactivation has no effect when
   * the process is already active.
   */
  bool only_one_per_event_;

  /**
   * Bias the dark brem cross section GLOBALLY
   */
  double global_bias_;

  /**
   * The mass of the A' during this run [MeV]
   */
  double ap_mass_;

  /**
   * Should we have a cache for the computed cross sections?
   */
  bool cache_xsec_;

  /**
   * Are we dark-brem off muons or electrons?
   */
  bool muons_;

  /**
   * The model that we are using in this run.
   *
   * Shared with the chaching class.
   */
  std::shared_ptr<G4DarkBremsstrahlungModel> model_;

  /// Our instance of a cross section cache
  ElementXsecCache element_xsec_cache_;

  /// the created process
  static G4DarkBremsstrahlung* the_process_;
};  // G4DarkBremsstrahlung

}  // namespace darkbrem
}  // namespace simcore

#endif  // SIMCORE_DARKBREM_G4DARKBREMSSTRAHLUNG_H_

