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

namespace parser {
class LHE {
  std::ifstream input_text_;
 public:
  LHE(const std::string& path)
    : input_text_{path} {
      if (not input_text_.open()) {
        throw std::runtime_error("Unable to open '"+path+"'.");
      }
  }
  ~LHE() {
    input_text_.close();
  }
  /**
   * Parse an LHE file from the input stream
   *
   * @param[in,out] lib library being constructed
   */
  void load(std::map<double, std::vector<OutgoingKinematics>>& lib) {
    static const double MA =
        G4APrime::APrime()->GetPDGMass() / CLHEP::GeV;  // mass A' in GeV
  
    std::string line;
    while (std::getline(is, line)) {
      std::istringstream iss(line);
      int ptype, state;
      double skip, px, py, pz, E, M;
      if (iss >> ptype >> state >> skip >> skip >> skip >> skip >> px >> py >>
          pz >> E >> M) {
        if ((ptype == 11 or ptype == 13) && (state == -1)) {
          double ebeam = E;
          double e_px, e_py, e_pz, a_px, a_py, a_pz, e_E, a_E, e_M, a_M;
          for (int i = 0; i < 2; i++) {
            std::getline(ifile, line);
          }
          std::istringstream jss(line);
          jss >> ptype >> state >> skip >> skip >> skip >> skip >> e_px >> e_py >>
              e_pz >> e_E >> e_M;
          if ((ptype == 11 or ptype == 13) && (state == 1)) {  // Find a final state lepton
            for (int i = 0; i < 2; i++) {
              std::getline(ifile, line);
            }
            std::istringstream kss(line);
            kss >> ptype >> state >> skip >> skip >> skip >> skip >> a_px >>
                a_py >> a_pz >> a_E >> a_M;
            if (ptype == 622 and state == 1) {
              if (abs(1. - a_M / MA) > 1e-3) {
                throw std::runtime_error(
                                "A MadGraph imported event has a different "
                                "APrime mass than the model has (MadGraph = " +
                                    std::to_string(a_M) + "GeV; Model = " +
                                    std::to_string(MA) + "GeV).");
              }
              OutgoingKinematics evnt;
              double cmpx = a_px + e_px;
              double cmpy = a_py + e_py;
              double cmpz = a_pz + e_pz;
              double cmE = a_E + e_E;
              evnt.lepton = LorentzVector(e_px, e_py, e_pz, e_E);
              evnt.centerMomentum = LorentzVector(cmpx, cmpy, cmpz, cmE);
              evnt.E = ebeam;
              lib[ebeam].push_back(evnt);
            }  // get a prime kinematics
          }    // check for final state
        }      // check for particle type and state
      }        // able to get momentum/energy numbers
    }          // while getting lines
  }
};  // LHE

class CSV {
  std::ifstream input_text_;
 public:
  CSV(const std::string& path)
    : input_text_{path} {
      if (not input_text_.open()) {
        throw std::runtime_error("Unable to open '"+path+"'.");
      }
  }
  ~CSV() {
    input_text_.close();
  }
  void load(std::map<double, std::vector<OutgoingKinematics>>& lib) {
  }
};

class LHEGZ {
  FILE* source_{NULL};
  z_stream strm_;
 public:
  LHEGZ(const std::string& path)
    : source_{fopen(path.c_str())} {
      if (source_ == NULL) {
        throw std::runtime_error("Unable to open '"+path+"'.");
      }
    }
  ~LHEGZ() {
    if (source_ != NULL) fclose(source_);
  }
  void load(std::map<double, std::vector<OutgoingKinematics>>& lib) {
  }
};

}  // namspace parser

}

void parseLibrary(const std::string& path, std::map<double, std::vector<OutgoingKinematics>>& lib) {
  // Assumptions:
  //  - Directory passed is a flat directory (no sub directories) containing LHE
  //  files
  //  - LHE files are events generated with the correct mass point
  //
  // A future improvement could be parsing a directory and only selecting LHE files
  // the contain dark photons of the configure mass. This has not been implemented
  // because there has been no reason to merge dark brem event libraries corresponding
  // to different mass points.

  if (hasEnding(path, ".gz")) {
    // zlib compressed data file
  } if (hasEnding(path, ".csv")) {
    ParseCSV(path, lib);  
  } else if (hasEnding(path, ".lhe")) {
    std::ifstream ifile;
    ifile.open(fname.c_str());
    if (!ifile) {
      throw std::runtime_error("Unable to open LHE file '"+fname+"'.");
    }
    ParseLHE(ifile, lib);
    ifile.close();
  } else {
    // assume directory of files
    DIR *dir;            // handle to opened directory
    struct dirent *ent;  // handle to entry inside directory
    if ((dir = opendir(path.c_str())) != NULL) {
      // directory can be opened
      while ((ent = readdir(dir)) != NULL) {
        std::string fp = path + '/' + std::string(ent->d_name);
        if (hasEnding(fp,".lhe") or hasEnding(fp, ".lhe.gz")) {
          // file ends in '.lhe' or '.lhe.gz'
          parseLibrary(fp, lib);
        }
      }
      closedir(dir);
    }
  }
}

}
