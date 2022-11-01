/**
 * @file G4APrime.cxx
 * @brief Class creating the A' particle in Geant.
 * @author Michael Revering, University of Minnesota
 */

#include "G4DarkBreM/G4APrime.h"

#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "globals.hh"

G4APrime* G4APrime::theAPrime = 0;

G4APrime* G4APrime::APrime() {
  if (!theAPrime) {
    throw std::runtime_error("Attempting to access the APrime particle before it has been initialized.");
  }

  return theAPrime;
}

void G4APrime::Initialize(double mass, int id) {
  if (theAPrime)
    throw std::runtime_error("Attempting to initialize the APrime particle more than once.");

  /**
   * Here are the properties of the formal Geant4 dark photon we define.
   *
   * Property | Value 
   * ---|---
   * short name | A^1 
   * mass | **configured**
   * mass width | 0 
   * electric charge | 0 
   * spin | 0
   * parity | 0
   * conjugation | 0
   * isospin | 0
   * isospin3 | 0
   * Gparity | 0
   * long name | APrime
   * lepton number | 0
   * baryon number | 0
   * PDG ID encoding | **configured**
   * is stable (no decay) | true
   * lifetime | -1 (i.e. no decay)
   * decay table | nullptr (i.e. no decay)
   */
  theAPrime = new G4APrime(
      "A^1" /* short name */, mass * MeV,
      0. /* mass width */, 0. /*electric charge */, 
      0 /* spin */, 0 /* parity */, 0 /* conjugation */, 
      0 /* isospine */, 0 /* isospin3 */, 0 /* G parity */, 
      "APrime" /* long name */, 0 /* lepton number */, 
      0 /* baryon number */, id, true /* stable? */, 
      -1 /* lifetime */, nullptr /* decay table */
      );
  return;
}

