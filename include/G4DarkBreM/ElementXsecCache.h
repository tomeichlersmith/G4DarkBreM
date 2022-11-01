#ifndef G4DarkBREM_ELEMENTXSECCACHE_H
#define G4DarkBREM_ELEMENTXSECCACHE_H

#include <memory>

#include "G4DarkBreM/PrototypeModel.h"

namespace g4db {

/**
 * The cache of already computed cross sections
 *
 * We make a specific class for the cache in order
 * to keep the key encoding/decoding process in a central
 * location.
 */
class ElementXsecCache {
 public:
  /**
   * Default constructor
   *
   * Does nothing interesting, but no model for calculating cross section has
   * been set.
   */
  ElementXsecCache() = default;

  /**
   * Constructor with a model to calculate the cross section.
   */
  ElementXsecCache(std::shared_ptr<PrototypeModel> model)
      : model_{model} {}

  /**
   * Get the value of the cross section for the input variables
   * and calculate the cross section if it wasn't calculated before.
   *
   * @throws std::runtime_error if no model is available for calculating cross sections
   * @param[in] energy Energy of incident lepton [MeV]
   * @param[in] A atomic mass of element [atomic mass units]
   * @param[in] Z atomic number of element [num protons]
   * @returns cross section corresponding to the input parameters (including
   * units Geant4 style)
   */
  G4double get(G4double energy, G4double A, G4double Z);

  /**
   * Stream the entire table into the output stream.
   *
   * @param[in,out] o ostream to write to
   */
  void stream(std::ostream& o) const;

  /**
   * Overload the streaming operator for ease
   *
   * @param[in] o ostream to write to
   * @param[in] c cache to write out
   * @returns modified ostream
   */
  friend std::ostream& operator<<(std::ostream& o, const ElementXsecCache c) {
    c.stream(o);
    return o;
  }

 private:
  /// The type for the key we use in the cache
  typedef unsigned long int key_t;

  /// The maximum value of A
  static const key_t MAX_A{1000};

  /// The maximum value for energy [MeV]
  static const key_t MAX_E{1500000};

  /**
   * Compute a key for the cache map
   * Generating a unique key _after_ making the energy [MeV] an integer.
   * The atomic mass (A) and charge (Z) are given by Geant4 as doubles as well,
   * so I cast them to integers before computing the key.
   *
   * This is what you would edit if you want a more/less find-grained cache
   * of Xsecs. Right now, since the internal unit of energy in Geant4 is MeV,
   * the cache is binned at the 1MeV scale.
   *
   * @param[in] energy Energy of incident lepton [MeV]
   * @param[in] A atomic mass of element [atomic mass units]
   * @param[in] Z atomic number of element [num protons]
   * @returns unsigned integer cache key for these three inputs
   */
  key_t computeKey(G4double energy, G4double A, G4double Z) const;

 private:
  /// the actual map from cache keys to calculated cross sections
  std::map<key_t, G4double> the_cache_;

  /// shared pointer to the model for calculating cross sections
  std::shared_ptr<PrototypeModel> model_;

};  // ElementXsecCache

}

#endif
