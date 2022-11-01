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

G4APrime* G4APrime::APrime(G4double theMass) {
  if (!theAPrime) {
    if (theMass < 0)
      throw std::runtime_error("Attempting to access the APrime particle before it has been initialized.");

    /**
     * Here are the properties of the formal Geant4 dark photon we define.
     *
     * Property | Value 
     * ---|---
     * short name | A^1 
     * mass | configured 
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
     * PDG ID encoding | 622
     * is stable (no decay) | true
     * lifetime | -1 (i.e. no decay)
     * decay table | nullptr (i.e. no decay)
     */
    const G4String& name = "A^1";
    G4double mass = theMass;
    G4double width = 0.;
    G4double charge = 0;
    G4int iSpin = 0;
    G4int iParity = 0;
    G4int iConjugation = 0;
    G4int iIsospin = 0;
    G4int iIsospin3 = 0;
    G4int gParity = 0;
    const G4String& pType = "APrime";
    G4int lepton = 0;
    G4int baryon = 0;
    G4int encoding = 622;    // PDG ID
    G4bool stable = true;    // stable - no decay
    G4double lifetime = -1;  // stable - no decay
    G4DecayTable* decaytable = 0;

    theAPrime =
        new G4APrime(name, mass, width, charge, iSpin, iParity, iConjugation,
                     iIsospin, iIsospin3, gParity, pType, lepton, baryon,
                     encoding, stable, lifetime, decaytable);
  }

  return theAPrime;
}

