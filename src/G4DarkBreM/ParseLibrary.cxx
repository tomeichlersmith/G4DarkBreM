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
 * We go line-by-line through the input text stream, looking for dark brem events.
 * A "dark brem event" in this context is defined below.
 *
 * ```
 *   lepton_id -1 <skip> <skip> <skip> <skip> px py pz E m
 *   <skip-line>
 *   <skip-line>
 *   lepton_id 1 <skip> <skip> <skip> <skip> px py pz E m
 *   <skip-line>
 *   <skip-line>
 *   aprime_id 1 <skip> <skip> <skip> <skip> px py pz E m
 * ```
 *
 * This matches a subcomponent of the LHE scheme written by MadGraph/MadEvent
 * (hence the reason this is the "lhe" parser); however, a lot of information
 * is skipped and additional assumptions are made in order to increase the 
 * parsing speed.
 *
 * The `lepton_id` is allowed to be _either_ 11 or 13 _everywhere_. No consistency
 * checking is done.
 *
 * The `E` from the first line is used as the incident lepton energy.
 * The four-momentum from the middle line is the recoil lepton's four momentum,
 * and the four-momentum from the last line is used in conjuction with the
 * recoil four-momentum to calculate the center of momentum vector.
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
    if (iss >> ptype >> state >> skip >> skip >> skip >> skip 
            >> px >> py >> pz >> E >> M) {
      if ((ptype == 11 or ptype == 13) && (state == -1)) {
        double incident_energy = E;
        double e_px, e_py, e_pz, a_px, a_py, a_pz, e_E, a_E, e_M, a_M;
        for (int i = 0; i < 2; i++) {
          std::getline(reader,line);
        }
        std::istringstream jss(line);
        jss >> ptype >> state >> skip >> skip >> skip >> skip 
            >> e_px >> e_py >> e_pz >> e_E >> e_M;
        if ((ptype == 11 or ptype == 13) && (state == 1)) {  // Find a final state lepton
          for (int i = 0; i < 2; i++) {
            std::getline(reader,line);
          }
          std::istringstream kss(line);
          kss >> ptype >> state >> skip >> skip >> skip >> skip 
              >> a_px >> a_py >> a_pz >> a_E >> a_M;
          if (ptype == aprime_lhe_id and state == 1) {
            OutgoingKinematics evnt;
            double cmpx = a_px + e_px;
            double cmpy = a_py + e_py;
            double cmpz = a_pz + e_pz;
            double cmE = a_E + e_E;
            evnt.lepton = CLHEP::HepLorentzVector(e_px, e_py, e_pz, e_E);
            evnt.centerMomentum = CLHEP::HepLorentzVector(cmpx, cmpy, cmpz, cmE);
            evnt.E = incident_energy;
            lib[incident_energy].push_back(evnt);
          }  // get a prime kinematics
        }    // check for final state
      }      // check for particle type and state
    }        // able to get momentum/energy numbers
  }          // while getting lines
}

/**
 * parse the input stream as a CSV file, filling the input library
 *
 * The CSV file is expected to have a **single** header line which
 * names the columns. These column names have no requirements
 * (besides the existence of this line).
 *
 * The CSV is required to have 9 columns on all non-empty lines of the file.
 * The 9 columns of the CSV all are in MeV and _in order_ are
 * 1. The incident lepton energy
 * 2. The total energy of the recoil
 * 3. The x-component of the recoil momentum
 * 4. The y-component of the recoil momentum
 * 5. The z-component of the recoil momentum
 * 6. The total energy of the A'
 * 7. The x-component of the A' momentum
 * 8. The y-component of the A' momentum
 * 9. The z-component of the A' momentum
 *
 * @note If developing this function, make sure to update dumpLibrary
 * so that they can be used in conjuction.
 *
 * @param[in] reader input stream reading the file
 * @param[in,out] lib dark brem event library to fill
 */
void csv(boost::iostreams::filtering_istream& reader, std::map<double, std::vector<OutgoingKinematics>>& lib) {
  std::string line;
  // skip the header line
  if (not std::getline(reader, line)) {
    throw std::runtime_error("Empty CSV file.");
  }
  // read in all non-empty lines
  while (std::getline(reader, line) and not line.empty()) {
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
  if (hasEnding(path, ".csv") or hasEnding(path, ".csv.gz") or hasEnding(path, ".lhe") or hasEnding(path, ".lhe.gz")) {
    /**
     * If the input path has one of the four file extensions below,
     * we assume it is a file to be parsed into the library.
     * - '.csv'
     * - '.csv.gz'
     * - '.lhe'
     * - '.lhe.gz'
     *
     * If the extension ends with '.gz', then a decompression step is
     * added to the input stream. Boost.Iostream provides the 
     * gzip_decompressor "filter" which does the decompression on the
     * data stream as it is being read in.
     *
     * @see parse::csv for files ending with '.csv' or '.csv.gz'
     * @see parse::lhe for files ending with '.lhe' or '.lhe.gz'
     */
    boost::iostreams::filtering_istream reader;
    if (hasEnding(path, ".gz")) reader.push(boost::iostreams::gzip_decompressor());
    reader.push(boost::iostreams::file_source(path)); 
    if (hasEnding(path, ".csv") or hasEnding(path, ".csv.gz")) parse::csv(reader, lib); 
    else parse::lhe(reader, aprime_lhe_id, lib);
  } else {
    /**
     * If the input path _does not_ match one of the four accepted
     * extensions, then we assume it is a directory 
     *
     * @note We _do not_ recursively enter subdirectories.
     */
    DIR *dir;            // handle to opened directory
    struct dirent *ent;  // handle to entry inside directory
    if ((dir = opendir(path.c_str())) != NULL) {
      // directory can be opened
      while ((ent = readdir(dir)) != NULL) {
        std::string fp = path + '/' + std::string(ent->d_name);
        if (hasEnding(fp,".lhe") or hasEnding(fp, ".lhe.gz")
            or hasEnding(fp,".csv") or hasEnding(fp, ".csv.gz")) {
          /**
           * If any of the directory entries has one of the acceptable
           * extensions, we recursively call this function on that 
           * file path so that it can be parsed into the library.
           */
          parseLibrary(fp, aprime_lhe_id, lib);
        }
      }
      closedir(dir);
    } else {
      /**
       * If we can't open the path that we assumed was a directory as a directory,
       * we end processing.
       */
      throw std::runtime_error("Unable to open '"+path+"' as a directory.");
    }
  }
}

void dumpLibrary(std::ostream& o, const std::map<double, std::vector<OutgoingKinematics>>& lib) {
  /**
   * This function writes out the input library as CSV to the input output stream
   * in the same format as expected by parse::csv.
   */
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
