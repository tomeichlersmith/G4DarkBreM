#ifndef G4DARKBREM_PARSELIBRARY_H
#define G4DARKBREM_PARSELIBRARY_H

namespace g4db {

/**
 * Data frame to store necessary information from LHE files
 */
struct OutgoingKinematics {
  /// 4-momentum of lepton in center of momentum frame for electron-A'
  /// system
  LorentzVector lepton;
  /// 4-vector pointing to center of momentum frame
  LorentzVector centerMomentum;
  /// energy of lepton before brem (used as key in mad graph data map)
  G4double E;
};  // OutgoingKinematics

/**
 * parse the input library and return the in-memory kinematics library
 *
 * @param[in] path path to library to parse
 * @param[in,out] lib map of incident energy keys to set of outgoing kinematics
 */
void parseLibrary(const std::string& path, std::map<double, std::vector<OutgoingKinematics>>& lib);

}

#endif
