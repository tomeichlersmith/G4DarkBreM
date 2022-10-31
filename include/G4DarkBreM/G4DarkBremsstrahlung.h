/**
 * @file G4DarkBremsstrahlung.h
 * @brief Class providing the Dark Bremsstrahlung process class.
 * @author Michael Revering, University of Minnesota
 * @author Tom Eichlersmith, University of Minnesota
 */

#ifndef G4DARKBREM_G4DARKBREMSSTRAHLUNG_H_
#define G4DARKBREM_G4DARKBREMSSTRAHLUNG_H_

// Geant
#include "G4VDiscreteProcess.hh"

#include "G4DarkBreM/PrototypeModel.h"
#include "G4DarkBreM/ElementXsecCache.h"

class G4String;
class G4ParticleDefinition;

/**
 * @class G4DarkBremsstrahlung
 *
 * Class that represents the dark brem process.
 * A muon or electron is allowed to brem a dark photon
 */
class G4DarkBremsstrahlung : public G4VDiscreteProcess {
 public:
  /**
   * The name of this process in Geant4
   */
  static const std::string PROCESS_NAME;

  /**
   * Constructor
   *
   * Configures this process by doing three main things:
   *  1. Registers this process with Geant4 as a 'fElectromagnetic' process
   *     - Needed for Geant4 Biasing Framework to recognize this process as
   * "bias-able"
   *  2. Defines the EM subtype as one different from all other EM processes
   *     - Needed so we don't replace another EM process
   *  3. Configures the process
   *
   * We add ourselves to the process table for muons (or electrons),
   * so the caller **does not** need to do this. Normally, this is done
   * in a physics constructor, but we have chosen to do this differently
   * to avoid mistakenly adding the dark brem process configured for muons
   * to the electron table (and vis-versa).
   *
   * @note Make sure to call SetModel immediately after constructing
   * the process.
   *
   * @param[in] muons true if muons are dark bremming, false for electrons
   * @param[in] only_one_per_event true if de-activating process after first dark brem
   * @param[in] global_bias bias xsec globally by this factor
   * @param[in] cache_xsec true if we should cache xsecs at the MeV level of precision
   */
  G4DarkBremsstrahlung(bool muons, bool only_one_per_event = false, 
      double global_bias = 1., bool cache_xsec = true);

  /**
   * Decide on the model for dark bremming
   *
   * This can only be called once per run and must be called **immediately**
   * after construction so that the process is configured properly.
   *
   * @param[in] the_model model to use for dark brem simulation
   */
  void SetModel(std::shared_ptr<g4db::PrototypeModel> the_model);

  /**
   * Destructor
   */
  virtual ~G4DarkBremsstrahlung() = default;

  /**
   * Checks if the passed particle should be able to do this process
   *
   * @return true if particle is the configured lepton
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
   * Get a reference to the cross section cache.
   *
   * Again, this method is public only to be available to the executable
   * that generates a cross section table and testing.
   * Do not use this unless you really know what you are doing.
   */
  g4db::ElementXsecCache& getCache() { return element_xsec_cache_; }

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
   * in the constructor.
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
   * @note If a user wants to use this option, the dark brem process
   * should be _always_ re-activated in the RunManager::TerminateOneEvent 
   * method. This reactivation has no effect when the process is already active.
   */
  bool only_one_per_event_;

  /**
   * Bias the dark brem cross section GLOBALLY
   */
  double global_bias_;

  /**
   * Should we have a cache for the computed cross sections?
   */
  bool cache_xsec_;

  /**
   * Are we dark-bremming off muons or electrons?
   */
  bool muons_;

  /**
   * The model that we are using in this run.
   *
   * Shared with the chaching class.
   */
  std::shared_ptr<g4db::PrototypeModel> model_;

  /// Our instance of a cross section cache
  g4db::ElementXsecCache element_xsec_cache_;
};  // G4DarkBremsstrahlung

#endif

