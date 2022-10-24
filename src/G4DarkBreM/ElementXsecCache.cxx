#include "G4DarkBreM/ElementXsecCache.h"

namespace g4db {

G4double ElementXsecCache::get(G4double energy, G4double A, G4double Z) {
  key_t key = computeKey(energy, A, Z);
  if (the_cache_.find(key) == the_cache_.end()) {
    if (model_.get() == nullptr) {
      throw std::runtime_error(
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

}
