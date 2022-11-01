#include <dirent.h>
#include <fstream>
#include <sstream>

#include <boost/iostreams/operations.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/close.hpp>

#include "G4DarkBreM/ParseLibrary.h"

namespace g4db {

/**
 * Check if the input string has the input ending
 *
 * @param[in] str string to test
 * @param[in] end ending to have
 * @return true if the last characters of str are end
 */
bool hasEnding(const std::string& str, const std::string& end) {
  if (str.length() >= end.length()) {
    return (str.compare(str.length() - end.length(), end.length(), end) == 0);
  } else {
    return false;
  }
}

/**
 * namespace holding implementation of library parsing
 */
namespace parse {

/**
 * Parse an LHE file from the input stream
 *
 * @param[in] reader input stream reading the file
 * @param[in] aprime_lhe_id ID number of the dark photon within the LHE file
 * @param[in,out] lib dark brem event library to fill
 */
void lhe(boost::iostreams::filtering_istream& reader, int aprime_lhe_id, std::map<double, std::vector<OutgoingKinematics>>& lib) {
  std::string line;
  while (std::getline(reader,line)) {
    std::istringstream iss(line);
    int ptype, state;
    double skip, px, py, pz, E, M;
    if (iss >> ptype >> state >> skip >> skip >> skip >> skip >> px >> py >>
        pz >> E >> M) {
      if ((ptype == 11 or ptype == 13) && (state == -1)) {
        double ebeam = E;
        double e_px, e_py, e_pz, a_px, a_py, a_pz, e_E, a_E, e_M, a_M;
        for (int i = 0; i < 2; i++) {
          std::getline(reader,line);
        }
        std::istringstream jss(line);
        jss >> ptype >> state >> skip >> skip >> skip >> skip >> e_px >> e_py >>
            e_pz >> e_E >> e_M;
        if ((ptype == 11 or ptype == 13) && (state == 1)) {  // Find a final state lepton
          for (int i = 0; i < 2; i++) {
            std::getline(reader,line);
          }
          std::istringstream kss(line);
          kss >> ptype >> state >> skip >> skip >> skip >> skip >> a_px >>
              a_py >> a_pz >> a_E >> a_M;
          if (ptype == aprime_lhe_id and state == 1) {
            OutgoingKinematics evnt;
            double cmpx = a_px + e_px;
            double cmpy = a_py + e_py;
            double cmpz = a_pz + e_pz;
            double cmE = a_E + e_E;
            evnt.lepton = CLHEP::HepLorentzVector(e_px, e_py, e_pz, e_E);
            evnt.centerMomentum = CLHEP::HepLorentzVector(cmpx, cmpy, cmpz, cmE);
            evnt.E = ebeam;
            lib[ebeam].push_back(evnt);
          }  // get a prime kinematics
        }    // check for final state
      }      // check for particle type and state
    }        // able to get momentum/energy numbers
  }          // while getting lines
}

/**
 * parse the input stream as a CSV file, filling the input library
 *
 * @param[in] reader input stream reading the file
 * @param[in,out] lib dark brem event library to fill
 */
void csv(boost::iostreams::filtering_istream& reader, std::map<double, std::vector<OutgoingKinematics>>& lib) {
  std::string line;
  if (not std::getline(reader_, line)) {
    throw std::runtime_error("Empty CSV file.");
  }
  while (std::getline(reader_, line)) {
    std::istringstream lss{line};
    std::vector<double> vals;
    std::string cell;
    while (std::getline(lss,cell,',')) {
      vals.push_back(std::stod(cell));
    }
    if (not lss and cell.empty()) vals.push_back(-9999);
    if (vals.size() != 9) {
      throw std::runtime_error("Malformed row in CSV file: not exactly 9 columns");
    }
    OutgoingKinematics ok;
    ok.E = vals[0];
    ok.lepton = CLHEP::HepLorentzVector(vals[2], vals[3], vals[4], vals[1]);
    ok.centerMomentum = CLHEP::HepLorentzVector(vals[6], vals[7], vals[8], vals[5]);
    lib[ok.E].push_back(ok);
  }
}

}  // namspace parser

void parseLibrary(const std::string& path, int aprime_lhe_id, std::map<double, std::vector<OutgoingKinematics>>& lib) {
  // Assumptions:
  //  - Directory passed is a flat directory (no sub directories) containing LHE
  //  files
  //  - LHE files are events generated with the correct mass point
  //
  // A future improvement could be parsing a directory and only selecting LHE files
  // the contain dark photons of the configure mass. This has not been implemented
  // because there has been no reason to merge dark brem event libraries corresponding
  // to different mass points.

  if (hasEnding(path, ".csv") or hasEnding(path, ".csv.gz") or hasEnding(path, ".lhe") or hasEnding(path, ".lhe.gz")) {
    // accepted file extensions
    boost::iostreams::filtering_istream reader;
    if (hasEnding(path, ".gz")) reader.push(boost::iostreams::gzip_decompressor());
    reader.push(boost::iostreams::file_source(path)); 
    if (hasEnding(path, ".csv") or hasEnding(path, ".csv.gz")) parse::csv(reader, lib); 
    else parse::lhe(reader, aprime_lhe_id, lib);
  } else {
    // assume directory of files
    DIR *dir;            // handle to opened directory
    struct dirent *ent;  // handle to entry inside directory
    if ((dir = opendir(path.c_str())) != NULL) {
      // directory can be opened
      while ((ent = readdir(dir)) != NULL) {
        std::string fp = path + '/' + std::string(ent->d_name);
        if (hasEnding(fp,".lhe") or hasEnding(fp, ".lhe.gz")
            or hasEnding(fp,".csv") or hasEnding(fp, ".csv.gz")) {
          // file ends in one of the acceptable endings
          parseLibrary(fp, aprime_lhe_id, lib);
        }
      }
      closedir(dir);
    } else {
      throw std::runtime_error("Unable to open '"+path+"' as a directory.");
    }
  }
}

void dumpLibrary(std::ostream& o, const std::map<double, std::vector<OutgoingKinematics>>& lib) {
  o << "incident_energy,recoil_energy,recoil_px,recoil_py,recoil_pz,"
         "centerMomentum_energy,centerMomentum_px,centerMomentum_py,centerMomentum_pz\n";
  for (const auto& lib_entry : lib) {
    for (const auto& sample : lib_entry.second) {
      o << sample.E << ','
        << sample.lepton.e() << ','
        << sample.lepton.px() << ','
        << sample.lepton.py() << ','
        << sample.lepton.pz() << ','
        << sample.centerMomentum.e() << ','
        << sample.centerMomentum.px() << ','
        << sample.centerMomentum.py() << ','
        << sample.centerMomentum.pz() << '\n';
    }
  }
  o.flush();
}

}
