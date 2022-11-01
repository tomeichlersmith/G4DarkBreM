#ifndef G4DARKBREM_G4DARKBREMMODEL_H
#define G4DARKBREM_G4DARKBREMMODEL_H

#include <memory>
#include <map>

#include "G4DarkBreM/PrototypeModel.h"

#include "CLHEP/Vector/LorentzVector.h"

namespace g4db {

/**
 * @class G4DarkBreMModel
 *
 * Geant4 implementation of the model for a particle undergoing a dark brem
 * where we use an imported event library to decide the outgoing kinematics.
 *
 * This is where all the heavy lifting in terms of calculating cross sections
 * and actually having a lepton do a dark brem occurs. This model depends
 * on several configurable parameters.
 *
 * - library_path : the full path to the directory containing the LHE dark brem
 *   vertices that will be read in to make the vertex library
 * - epsilon : strength of the dark photon - photon mixing
 * - threshold : minimum energy in GeV for the lepton to have a non-zero
 *   cross section for going dark brem
 * - method : scaling method to use to scale the dark brem vertices from
 *   the library to the actual lepton energy when a dark brem occurs
 * - whether we are dark bremming of muons or electrons
 *
 * The required parameter is a vertex library generated in MadGraph
 * (library_path).
 */
class G4DarkBreMModel : public PrototypeModel {
  /// the class we use to store four vectors
  using LorentzVector = CLHEP::HepLorentzVector;
 public:
  /**
   * Set the parameters for this model.
   *
   * @param[in] method_name converted to an enum through a hard-coded switch statement.
   * @param[in] threshold minimum energy lepton needs to have to dark brem [GeV]
   * @param[in] epsilon dark photon mixing strength
   * @param[in] library_path directory in which MG library is stored
   * @param[in] muons true if using muons, false for electrons
   * @param[in] load_library only used in cross section executable where it is known
   *            that the library will not be used during program run
   *
   * The threshold is set to the maximum of the passed value or twice
   * the A' mass (so that it kinematically makes sense).
   *
   * The library path is immediately passed to SetMadGraphDataLibrary.
   */
  G4DarkBreMModel(const std::string& method_name, double threshold, 
      double epsilon, const std::string& library_path, bool muons, 
      bool load_library = true);

  /**
   * Destructor
   */
  virtual ~G4DarkBreMModel() = default;

  /**
   * Print the configuration of this model
   */
  virtual void PrintInfo() const;

  /**
   * Calculates the cross section per atom in GEANT4 internal units.
   *
   * The estimate for the total cross section given the material and the lepton's energy is done using an
   * implementation of the WW approximation using Boost's Math Quadrature library to numerically calculate 
   * the integrals. The actual formulas are listed here for reference.
   *
   * Since muons and electrons have such different masses, different approaches were required to create an 
   * approximation that both follows the trend produced by MG/ME and is sufficiently quick.
   *
   * ## Electrons
   * Because the electron mass is small, it typically suffices to calculate the effective photon flux \f$\chi\f$ 
   * once rather than modeling its functional dependence on the \f$\aprime\f$ energy and angle, as in the "full" 
   * WW approximation used below for muons. With electron's low mass, the Improved WW approximation can be used:
   * \f{equation}{
   * \sigma = \frac{pb}{GeV} \chi \int_0^{\min(1-m_e/E_0,1-m_A/E_0)} \frac{d\sigma}{dx}(x)dx
   * \f}
   * where
   * \f{equation}{
   * \chi = \int^{m_A^2}_{m_A^4/(4E_0^2)} dt \left( \frac{Z^2a^4t^2}{(1+a^2t)^2(1+t/d)^2}+\frac{Za_p^4t^2}{(1+a_p^2t)^2(1+t/0.71)^8}\left(1+\frac{t(\mu_p^2-1)}{4m_p^2}\right)^2\right)\frac{t-m_A^4/(4E_0^2)}{t^2}
   * \f}
   * \f{equation}{
   * a = \frac{111.0}{m_e Z^{1/3}}
   * \quad
   * a_p = \frac{773.0}{m_e Z^{2/3}}
   * \quad
   * d = \frac{0.164}{A^{2/3}}
   * \f}
   * \f{equation}{
   * \frac{d\sigma}{dx}(x) = 4 \alpha_{EW}^3\epsilon^2 \sqrt{1-\frac{m_A^2}{E_0^2}}\frac{1-x+x^2/3}{m_A^2(1-x)/x+m_e^2x}
   * \f}
   *
   * - \f$E_0\f$ is the incoming electrons's energy in GeV
   * - \f$m_e\f$ is the mass of the electron in GeV
   * - \f$m_A\f$ is the mass of the dark photon in GeV
   * - \f$m_p = 0.938\f$ is the mass of the proton in GeV
   * - \f$\mu_p = 2.79\f$ is the proton \f$\mu\f$
   * - \f$A\f$ is the atomic mass of the target nucleus in amu
   * - \f$Z\f$ is the atomic number of the target nucleus
   * - \f$\alpha_{EW} = 1/137\f$ is the fine-structure constant
   * - \f$\epsilon\f$ is the dark photon mixing strength the the SM photon
   * - \f$pb/GeV = 3.894\times10^8\f$ is a conversion factor from GeV to pico-barns.
   * 
   * ## Muons 
   * The muon's greater mass motivated the use of the "full" WW but including the numerical evaluation of 
   * \f$\chi\f$ at each point in phase space proved to be too costly. 
   * Instead, we use an analytic integration of only elastic form-factor component:
   * 
   * \f{equation}{
   * \sigma = \frac{pb}{GeV} \int_0^{0.3} \int_0^{\min(1-m_\mu/E_0,1-m_A/E_0)} \frac{d\sigma}{dxd\theta}~dx~d\theta
   * \f}
   * 
   * where
   * 
   * \f{equation}{
   * \frac{d\sigma}{dx~d\cos\theta} = 2 \alpha_{EW}^3\epsilon^2 \sqrt{x^2E_0^2 - m_A^2}E_0(1-x)
   *     \frac{\chi(x,\theta)}{\tilde{u}^2} \mathcal{A}^2
   * \f}
   *
   * and
   *
   * \f{equation}{
   * \chi(x,\theta) = - \frac{Z^2(a^{-2}+d+2t_{max})}{(a^{-2}-d)^3}\left(
   *     \frac{(a^{-2}-d)(t_{max}-t_{min})}{(a^{-2}+t_{max})(d+t_{max})}
   *     + \log\left(\frac{(a^{-2}+t_{max})(d+t_{min})}{(a^{-2}+t_{min})(d+t_{max})}\right)
   *     \right)
   * \f}
   * \f{equation}{
   * \mathcal{A}^2 = 2\frac{2-2x+x^2}{1-x}+\frac{4(m_A^2+2m_\mu^2)}{\tilde{u}^2}(\tilde{u}x + m_A^2(1-x) + m_\mu^2x^2)
   * \f}
   * \f{equation}{
   * \tilde{u} = -xE_0^2\theta^2 - m_A^2\frac{1-x}{x} - m_\mu^2x
   * \f}
   * \f{equation}{
   * t_{min} = \left(\frac{\tilde{u}}{2E_0(1-x)}\right)^2 \qquad t_{max} = E_0^2
   * \f}
   * and \f$m_\mu\f$ is the mass of the muon in GeV, and the other symbols are the same as the electron case.
   *
   * @param lepton_ke kinetic energy of incoming particle
   * @param atomicZ atomic number of atom
   * @param atomicA atomic mass of atom
   * @return cross section (0. if outside energy cuts)
   */
  virtual G4double ComputeCrossSectionPerAtom(G4double lepton_ke,
                                              G4double atomicA,
                                              G4double atomicZ);

  /**
   * Scale one of the MG events in our library to the input incident 
   * lepton energy.
   *
   * This is also helpful for testing the scaling procedure in its own
   * executable separate from the Geant4 infrastructure.
   *
   * @note The vector returned is relative to the incident lepton as if
   * it came in along the z-axis.
   *
   * Gets an energy fraction and transverse momentum (\f$p_T\f$) from the
   * loaded library of MadGraph events using the entry in the library with
   * the nearest incident energy above the actual input incident energy.
   *
   * The scaling of this energy fraction and \f$p_T\f$ to the actual lepton
   * energy depends on the input method. In all cases, the azimuthal angle
   * is chosen uniformly between 0 and \f$2\pi\f$.
   *
   * ## Forward Only
   * Scales the energy so that the fraction of kinetic energy is constant,
   * keeping the \f$p_T\f$ constant. 
   *
   * If the \f$p_T\f$ is larger than the new energy, that event
   * is skipped, and a new one is taken from the file. If the loaded library
   * does not fully represent the range of incident energies being seen
   * by the simulation, this will occur frequently.
   *
   * With only the kinetic energy fraction and \f$p_T\f$, the sign of
   * the longitudinal momentum \f$p_z\f$ is undetermined. This method
   * simply chooses the \f$p_z\f$ of the recoil lepton to always be positive.
   *
   * ## CM Scaling
   * Scale MadGraph vertex to actual energy of lepton using Lorentz boosts.
   *
   * The scaling is done via two boosts.
   * 1. Boost out of the center-of-momentum (CoM) frame read in along with the
   *    MadGraph event library.
   * 2. Boost into approximately) the incident lepton energy frame by
   *    constructing a "new" CoM frame using the actual CoM frame's 
   *    transverse momentum and lowering the \f$p_z\f$ and energy of 
   *    the CoM by the difference between the input incident
   *    energy and the sampled incident energy.
   *
   * After these boosts, the energy of the recoil and its \f$p_T\f$ are
   * extracted.
   *
   * ## Undefined
   * Don't scale the MadGraph vertex to the actual energy of the lepton.
   *
   * We simply copy the read-in recoil energy's energy, momentum, and \f$p_T\f$.
   *
   * @param[in] incident_energy incident total energy of the lepton [GeV] 
   * @param[in] lepton_mass mass of incident lepton [GeV]
   * @return G4ThreeVector representing the recoil lepton's outgoing momentum
   */
  G4ThreeVector scale(double incident_energy, double lepton_mass);

  /**
   * Simulates the emission of a dark photon + lepton
   *
   * @see scale for how the event library is sampled and scaled to the incident
   * lepton's actual energy
   *
   * After calling scale, we rotate the outgoing lepton's momentum to the
   * frame of the incident particle and then calculate the dark photon's
   * momentum such that three-momentum is conserved.
   *
   * @param[in,out] particleChange structure holding changes to make to particle
   * track
   * @param[in] track current track being processesed
   * @param[in] step current step of the track
   */
  virtual void GenerateChange(G4ParticleChange& particleChange,
                              const G4Track& track, const G4Step& step);

  /**
   * @struct OutgoingKinematics
   *
   * Data frame to store mad graph data read in from LHE files.
   */
  struct OutgoingKinematics {
    /// 4-momentum of lepton in center of momentum frame for electron-A'
    /// system
    LorentzVector lepton;
    /// 4-vector pointing to center of momentum frame
    LorentzVector centerMomentum;
    /// energy of lepton before brem (used as key in mad graph data map)
    G4double E;
  };
 private:
  /**
   * Set the library of dark brem events to be scaled.
   *
   * This function loads the directory of LHE files passed
   * into our in-memory library of events to be sampled from.
   *
   * @param path path to directory of LHE files
   */
  void SetMadGraphDataLibrary(std::string path);

  /**
   * Parse an LHE File
   *
   * Parses an LHE file to extract the kinetic energy fraction and pt of the
   * outgoing electron in each event. Loads the two numbers from every event
   * into a map of vectors of pairs (mgdata). Map is keyed by energy, vector
   * pairs are energy fraction + pt. Also creates an list of energies and
   * placeholders (energies), so that different energies can be looped
   * separately.
   *
   * @param fname name of LHE file to parse
   */
  void ParseLHE(std::string fname);

  /**
   * Fill vector of currentDataPoints_ with the same number of items as the
   * madgraph data.
   *
   * Randomly choose a starting point so that the simulation run isn't dependent
   * on the order of the events as written in the LHE library. The random starting
   * position is uniformly chosen using G4Uniform() so. The sample function
   * will loop from the last event parsed back to the first event parsed so that
   * the starting position does not matter.
   *
   * In this function, we also update maxIterations_ so that it is equal to the smallest
   * entry in the library (with a maximum of 10k). This saves time in the situation where
   * an incorrect library was accidentally used and the simulation is looping through events
   * attempting to find one that can fit its criteria.
   */
  void MakePlaceholders();

  /**
   * Returns MadGraph data given an energy [GeV].
   *
   * Gets the energy fraction and \f$p_T\f$ from the imported LHE data.
   * incident_energy should be in GeV, returns the sample outgoing kinematics.
   *
   * Samples from the closest imported incident energy _above_ the given value
   * (this helps avoid biasing issues).
   *
   * @param incident_energy energy of particle undergoing dark brem [GeV]
   * @return sample outgoing kinematics
   */
  OutgoingKinematics sample(double incident_energy);

 private:
  /**
   * maximum number of iterations to check before giving up on an event
   *
   * This is only used in the ForwardOnly scaling method and is only
   * reached if the event library energies are not appropriately matched
   * with the energy range of particles that are existing in the simulation.
   */
  unsigned int maxIterations_{10000};

  /** 
   * Threshold for non-zero xsec [GeV]
   *
   * Configurable with 'threshold'. At minimum, it is always
   * at least twice the dark photon mass.
   */
  double threshold_;

  /** 
   * Epsilon value to plug into xsec calculation
   *
   * @sa ComputeCrossSectionPerAtom for how this is used
   *
   * Configurable with 'epsilon'
   */
  double epsilon_;

  /**
   * @enum DarkBremMethod
   *
   * Possible methods to use the dark brem vertices from the imported library
   * inside of this model.
   */
  enum DarkBremMethod {
    /// Use actual lepton energy and get pT from LHE 
    /// (such that \f$p_T^2+m_l^2 < E_{acc}^2\f$)
    ForwardOnly = 1,
    /// Boost LHE vertex momenta to the actual lepton energy
    CMScaling = 2,
    /// Use LHE vertex as is
    Undefined = 3
  };

  /** method for this model
   *
   * Configurable with 'method'
   */
  DarkBremMethod method_{DarkBremMethod::Undefined};

  /**
   * Name of method for persisting into the RunHeader
   */
  std::string method_name_;

  /**
   * Full path to the vertex library used for persisting into the RunHeader
   */
  std::string library_path_;

  /**
   * should we always create a totally new lepton when we dark brem?
   *
   * @note make this configurable? I (Tom E) can't think of a reason NOT to have
   * it... The alternative is to allow Geant4 to decide when to make a new
   * particle by checking if the resulting kinetic energy is below some
   * threshold.
   */
  bool alwaysCreateNewLepton_{true};

  /**
   * Storage of data from mad graph
   *
   * Maps incoming lepton energy to various options for outgoing kinematics.
   * This is a hefty map and is what stores **all** of the events
   * imported from the LHE library of dark brem events.
   */
  std::map<double, std::vector<OutgoingKinematics> > madGraphData_;

  /**
   * Stores a map of current access points to mad graph data.
   *
   * Maps incoming lepton energy to the index of the data vector
   * that we will get the data from.
   *
   * Also sorts the incoming lepton energy so that we can find
   * the sampling energy that is closest above the actual incoming energy.
   */
  std::map<double, unsigned int> currentDataPoints_;
};

}  // namespace g4db

#endif  // G4DARKBREM_G4DARKBREMMODEL_H
