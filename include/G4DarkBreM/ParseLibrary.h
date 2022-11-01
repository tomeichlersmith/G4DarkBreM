/**
 * @file ParseLibrary.h
 * Declaration of library parsing function
 */

#ifndef G4DARKBREM_PARSELIBRARY_H
#define G4DARKBREM_PARSELIBRARY_H

#include <string>
#include <map>
#include <vector>
#include <ostream>

#include "CLHEP/Vector/LorentzVector.h"

namespace g4db {

/**
 * Data frame to store necessary information from LHE files
 */
struct OutgoingKinematics {
  /// 4-momentum of lepton in center of momentum frame for electron-A'
  /// system
  CLHEP::HepLorentzVector lepton;
  /// 4-vector pointing to center of momentum frame
  CLHEP::HepLorentzVector centerMomentum;
  /// energy of lepton before brem (used as key in mad graph data map)
  double E;
};  // OutgoingKinematics

/**
 * parse the input library and return the in-memory kinematics library
 *
 * @param[in] path path to library to parse
 * @param[in,out] lib map of incident energy keys to set of outgoing kinematics
 */
void parseLibrary(const std::string& path, int aprime_lhe_id, std::map<double, std::vector<OutgoingKinematics>>& lib);

/**
 * Dump the input library to the input output stream
 *
 * @note This is only helpful for our extraction executable
 * and for testing that the parsing is performing correctly.
 *
 * @param[in,out] o output stream to write CSV to
 * @param[in] lib library to write out
 */
void dumpLibrary(std::ostream& o, const std::map<double, std::vector<OutgoingKinematics>>& lib);

}

#endif
