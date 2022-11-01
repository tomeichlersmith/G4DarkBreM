/**
 * @file G4APrime.h
 * @brief Class creating the A' particle in Geant.
 * @author Michael Revering, University of Minnesota
 */

#ifndef SIMCORE_DARKBREM_G4APRIME_H_
#define SIMCORE_DARKBREM_G4APRIME_H_

// Geant
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"

class G4String;
class G4DecayTable;

/**
 * Formal class representing the A' (a dark photon)
 *
 * This class follows the standard prototype for all G4ParticleDefinitions
 * where a static private member holds onto the single instance of this
 * particle definition to be shared by everyone. In this case, 
 * _the first call_ to the accessor function needs to define the A' mass
 * and then all subsequent calls can access it like normal with that
 * configured mass.
 *
 * This first call to the accessing method should be done in the 
 * ConstructParticle function of a physics constructor.
 */
class G4APrime : public G4ParticleDefinition {
 private:
  /** Reference to single particle definition of A' */
  static G4APrime* theAPrime;

  /**
   * Constructor
   *
   * Passes all parameters to the base class constructor
   * to register this particle definition with Geant4.
   */
  G4APrime(const G4String& Name, G4double mass, G4double width, G4double charge,
           G4int iSpin, G4int iParity, G4int iConjugation, G4int iIsospin,
           G4int iIsospin3, G4int gParity, const G4String& pType, G4int lepton,
           G4int baryon, G4int encoding, G4bool stable, G4double lifetime,
           G4DecayTable* decaytable)
      : G4ParticleDefinition(Name, mass, width, charge, iSpin, iParity,
                             iConjugation, iIsospin, iIsospin3, gParity, pType,
                             lepton, baryon, encoding, stable, lifetime,
                             decaytable) {}

  /**
   * Destructor
   *
   * Does nothing on purpose.
   */
  virtual ~G4APrime() {}

 public:
  /**
   * Accessor for APrime definition
   *
   * @throws std::runtime_error if the APrime is not initialized yet.
   *
   * @see Initialize for configuring and constructing the APrime
   * at the start of a run.
   */
  static G4APrime* APrime();

  /**
   * Initialize the APrime particle with the passed configuration
   * 
   * @throws std::runtime_error if the APrime has already been initialized
   *
   * @param[in] mass The mass of the APrime in MeV
   * @param[in] id The PDG ID number to use for the APrime particle
   *
   * The default value for the PDG ID is set to 62 and has been arbitrarily
   * chosen from the range defined by 11(c) in the 
   * [PDG ID numbering scheme](https://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf),
   * avoiding the already-defined "one-of-a-kind" particles using 39, 41, and 42.
   */
  static void Initialize(double mass, int id = 62);
};

#endif  // SIMCORE_DARKBREM_G4APRIME_H_
