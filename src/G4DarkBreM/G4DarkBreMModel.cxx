
#include "G4DarkBreM/G4DarkBreMModel.h"
#include "G4DarkBreM/G4APrime.h"
#include "G4DarkBreM/ParseLibrary.h"

// Geant4
#include "Randomize.hh"
#include "G4Electron.hh"
#include "G4MuonMinus.hh"
#include "G4EventManager.hh"  //for EventID number
#include "G4PhysicalConstants.hh"
#include "G4RunManager.hh"  //for VerboseLevel
#include "G4SystemOfUnits.hh"

// Boost
#include <boost/math/quadrature/gauss_kronrod.hpp>

// STL
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

namespace g4db {

/**
 * The integration method we will be using for our numerical integrals
 *
 * The Gauss-Kronrod method was chosen due to its ability to limit the
 * number of calls to the function representing the integrand which
 * should help improve performance for us due to the complexity of our
 * integrand. The order of the GK method was chosen after some 
 * experimentation, starting at a high value (61) and then lowering
 * it to achieve better performance while checking the accuracy of
 * the results.
 *
 * As explained in the [Boost GK Docs](https://www.boost.org/doc/libs/master/libs/math/doc/html/math_toolkit/gauss_kronrod.html),
 * generally the error estimation technique for this method is
 * overly pessimistic, so we can confidently set the maximum
 * depth low and the desired relative error high compared
 * to other methods. We have followed the examples in the docs
 * where we use max_depth to 5 and relative error to 1e-9.
 */
using int_method = boost::math::quadrature::gauss_kronrod<double, 61>;

/**
 * numerically integrate the value of the flux factory chi
 *
 * The integration of the form factor into the flux factor can
 * be done analytically with a tool like mathematica, but when
 * including the inelastic term, it produces such a complicated 
 * result that the numerical integration is actually *faster*
 * than the analytical one.
 *
 * The form factors are copied from Appendix A (Eq A18 and A19) of
 * https://journals.aps.org/prd/pdf/10.1103/PhysRevD.80.075018
 */
static double flux_factor_chi_numerical(G4double A, G4double Z, double tmin, double tmax) {
  /*
   * bin = (mu_p^2 - 1)/(4 m_pr^2)
   * mel = mass of electron in GeV
   */
  static const double bin = (2.79*2.79 - 1)/(4*0.938*0.938),
                      mel = 0.000511;
  const double ael = 111.0*pow(Z,-1./3.)/mel,
               del = 0.164*pow(A,-2./3.),
               ain = 773.0*pow(Z,-2./3.)/mel,
               din = 0.71,
               ael_inv2 = pow(ael, -2),
               ain_inv2 = pow(ain, -2);

  /**
   * We've manually expanded the integrand to cancel out the 1/t^2 factor
   * from the differential, this helps the numerical integration converge
   * because we aren't teetering on the edge of division by zero
   *
   * The `auto` used in the integrand definition represents a _function_ 
   * whose return value is a `double` and which has a single input `t`. 
   * This lambda expression saves us the time of having to re-calculate 
   * the form factor constants that do not depend on `t` because it 
   * can inherit their values from the environment. 
   * The return value is a double since it is calculated
   * by simple arithmetic operations on doubles.
   */
  auto integrand = [&](double t) {
    double ael_factor = 1./(ael_inv2 + t),
           del_factor = 1./(1+t/del),
           ain_factor = 1./(ain_inv2 + t),
           din_factor = 1./(1+t/din),
           nucl = (1 + t*bin);
    
    return (pow(ael_factor*del_factor*Z, 2)
            +
            Z*pow(ain_factor*nucl*din_factor*din_factor*din_factor*din_factor, 2)
           )*(t-tmin);
  };

  return int_method::integrate(integrand,tmin,tmax,5,1e-9);
}

/**
 * analytic flux factor chi integrated and simplified by DMG4 authors
 *
 * This only includes the elastic form factor term
 */
static double flux_factor_chi_analytic(G4double A, G4double Z, double tmin, double tmax) {
  static const double mel = 0.000511;
  const double a_el = 111.*pow(Z,-1./3)/mel,
               d_el = 0.164*pow(A,-2./3);
  double ta = 1.0/(a_el*a_el);
  double td = d_el;
  return -Z*Z*((td*td*(
              ((ta - td)*(ta + td + 2.0*tmax)*(tmax - tmin))/((ta + tmax)*(td + tmax)) 
              + (ta + td + 2.0*tmin)*(log(ta + tmax) - log(td + tmax) - log(ta + tmin) + log(td + tmin))
             ))/((ta-td)*(ta-td)*(ta-td)));
}

G4DarkBreMModel::G4DarkBreMModel(const std::string& method_name, double threshold,
    double epsilon, const std::string& library_path, bool muons, int aprime_lhe_id, 
    bool load_library)
    : PrototypeModel(muons), maxIterations_{10000}, 
      threshold_{std::max(threshold, 2.*G4APrime::APrime()->GetPDGMass()/CLHEP::GeV)},
      epsilon_{epsilon}, aprime_lhe_id_{aprime_lhe_id}, 
      method_(DarkBremMethod::Undefined), method_name_{method_name}, 
      library_path_{library_path} {
  if (method_name_ == "forward_only") {
    method_ = DarkBremMethod::ForwardOnly;
  } else if (method_name_ == "cm_scaling") {
    method_ = DarkBremMethod::CMScaling;
  } else if (method_name_ == "undefined") {
    method_ = DarkBremMethod::Undefined;
  } else {
    throw std::runtime_error("Invalid dark brem interpretaion/scaling method '"+method_name_+"'.");
  }

  if (load_library) SetMadGraphDataLibrary(library_path_);
}

void G4DarkBreMModel::PrintInfo() const {
  G4cout << " Dark Brem Vertex Library Model" << G4endl;
  G4cout << "   Threshold [GeV]: " << threshold_ << G4endl;
  G4cout << "   Epsilon:         " << epsilon_ << G4endl;
  G4cout << "   Scaling Method:  " << method_name_ << G4endl;
  G4cout << "   Vertex Library:  " << library_path_ << G4endl;
}

G4double G4DarkBreMModel::ComputeCrossSectionPerAtom(
    G4double lepton_ke, G4double A, G4double Z) {
  static const double MA = G4APrime::APrime()->GetPDGMass() / GeV;
  static const double MA2 = MA*MA;
  static const double alphaEW = 1.0 / 137.0;

  const double lepton_mass{
    (muons_ ? G4MuonMinus::MuonMinus()->GetPDGMass() : G4Electron::Electron()->GetPDGMass()) / GeV};
  const double lepton_mass_sq{lepton_mass*lepton_mass};

  // the cross section is zero if the lepton does not have enough
  // energy to create an A'
  // the threshold_ can also be set by the user to a higher value
  // to prevent dark-brem within inaccessible regions of phase
  // space
  if (lepton_ke < keV or lepton_ke < threshold_*GeV) return 0.;

  // Change energy to GeV.
  double lepton_e = lepton_ke/GeV + lepton_mass;
  double lepton_e_sq = lepton_e*lepton_e;

  /*
   * "Hyper-Improved" WW
   *
   * assume theta = 0, and x = 1 for form factor integration
   * i.e. now chi is a constant pulled out of the integration
   */
  double chi_hiww = flux_factor_chi_numerical(A,Z,MA2*MA2/(4*lepton_e_sq),MA2+lepton_mass_sq);

  /*
   * Differential cross section with respect to x and theta
   *
   * Equation (16) from Appendix A of https://arxiv.org/pdf/2101.12192.pdf
   *
   * This `auto` represents a lambda-expression function, inheriting many
   * pre-calculated constants (like lepton_e and chi) while also calculating
   * the variables dependent on the integration variables. The return value
   * of this function is a double since it is calculated by arithmetic
   * operations on doubles.
   */
  auto diff_cross = [&](double x, double theta) {
    if (x*lepton_e < threshold_) return 0.;

    double theta_sq = theta*theta;
    double x_sq = x*x;

    double utilde = -x*lepton_e_sq*theta_sq - MA2*(1.-x)/x - lepton_mass_sq*x;
    double utilde_sq = utilde*utilde;

    /*
     * WW
     *
     * Since muons are so much more massive than electrons, we keep 
     * the form factor integration limits dependent on x and theta
     */

    // non-zero theta and non-zero m_l
    double tmin = utilde_sq/(4.0*lepton_e_sq*(1.0-x)*(1.0-x));
    // maximum t kinematically limited to the incident lepton energy
    double tmax = lepton_e_sq;

    /*
     * The chi integrand limits given by
     *
     * Eqs (3.20) and (A6) of
     * https://journals.aps.org/prd/pdf/10.1103/PhysRevD.8.3109
     * OR
     * Eqs (3.2) and (3.6) of 
     * https://journals.aps.org/rmp/pdf/10.1103/RevModPhys.46.815
     *
     * to be
     *
     * tmax = m^2(1+l)^2
     * tmin = tmax / (2*E*x*(1-x))^2
     *
     * where
     *
     *  l = E^2x^2theta^2/m^2
     *  m is mass of dark photon
     *  E is the incident lepton energy
     * 
     * were investigated in an attempt to control the numerical integration
     * of chi in the hopes that cutting the integral away from odd places
     * would be able to avoid the funky business. This was not successful,
     * but we are leaving them here in case a typo is found in the future
     * or the search is chosen to resume.
    double el = lepton_e_sq*x_sq*theta_sq/MA2;
    double tmax = MA2*pow(1 + el,2);
    double tmin = tmax / pow(2*lepton_e*x*(1-x),2);
     */
  
    // require 0 < tmin < tmax to procede
    if (tmin < 0) return 0.;
    if (tmax < tmin) return 0.;
  
    /*
     * numerically integrate to calculate chi ourselves
     * this _has not_ been well behaved due to the extreme values
     * of t that must be handled
    double chi = flux_factor_chi_numerical(A,Z, tmin, tmax);
     */
  
    /*
     * use analytic elastic-only chi derived for DMG4
     * and double-checked with Mathematica
     *
     * The inelastic integral contains some 4000 terms
     * according to Mathematica so it is expensive to
     * compute and only an O(few) percent change.
     */
    double chi_analytic_elastic_only = flux_factor_chi_analytic(A,Z,tmin,tmax);
    
    /*
     * Amplitude squared is taken from 
     * Equation (17) from Appendix A of https://arxiv.org/pdf/2101.12192.pdf
     * with X = V
     */
    double factor1 = 2.0*(2.0 - 2.*x + x_sq)/(1. - x);
    double factor2 = 4.0*(MA2 + 2.0*lepton_mass_sq)/utilde_sq;
    double factor3 = utilde*x + MA2*(1. - x) + lepton_mass_sq*x_sq;
    double amplitude_sq = factor1 + factor2*factor3;

    return 2.*pow(epsilon_,2.)*pow(alphaEW,3.)
             *sqrt(x_sq*lepton_e_sq - MA2)*lepton_e*(1.-x)
             *(chi_analytic_elastic_only/utilde_sq)*amplitude_sq*sin(theta);
  };

  // deduce integral bounds
  double xmin = 0;
  double xmax = 1;
  if ((lepton_mass / lepton_e) > (MA / lepton_e))
    xmax = 1 - lepton_mass / lepton_e;
  else
    xmax = 1 - MA / lepton_e;

  /*
   * max recoil angle of A'
   *
   * The wide angle A' are produced at a negligible rate
   * so we enforce a hard-coded cut-off to stay within
   * the small-angle regime.
   *
   * We choose the same cutoff as DMG4.
   */
  double theta_max{0.3};

  /*
   * Integrand for integral over x
   *
   * For muons, we want to include the variation over theta from the chi
   * integral, so we calculate the x-integrand by numerically integrating
   * over theta in the differential cross section defined above.
   *
   * For electrons, we are using the Improved WW method where the theta
   * integral has already been done analytically and we can use the
   * numerical Chi (including both inelastic and elastic form factors)
   * calculated above.  
   *
   * This is the final lambda expression used here. Its one argument is a double
   * and it returns a double.
   */
  auto theta_integral = [&](double x) {
    if (muons_) {
      auto theta_integrand = [&](double theta) {
        return diff_cross(x, theta);
      };
      // integrand, min, max, max_depth, tolerance, error, pL1
      return int_method::integrate(theta_integrand, 0., theta_max, 5, 1e-9);
    } else {
      if (x*lepton_e < threshold_) return 0.;
      double beta = sqrt(1 - MA2/lepton_e_sq),
             nume = 1. - x + x*x/3.,
             deno = MA2*(1-x)/x + lepton_mass_sq;
      return 4*pow(epsilon_,2)*pow(alphaEW,3)*chi_hiww*beta*nume/deno;
    }
  };

  double error;
  double integrated_xsec = int_method::integrate(theta_integral, xmin, xmax, 5, 1e-9, &error);

  G4double GeVtoPb = 3.894E08;

  /*
   * The integrated_xsec should be the correct value, we are just
   * converting it to Geant4's pb units here
   */
  G4double cross = integrated_xsec * GeVtoPb * CLHEP::picobarn;

  if (cross < 0.) return 0.;  // safety check all the math

  return cross;
}

G4ThreeVector G4DarkBreMModel::scale(double incident_energy, double lepton_mass) {
  // mass A' in GeV
  static const double MA = G4APrime::APrime()->GetPDGMass() / CLHEP::GeV;
  OutgoingKinematics data = sample(incident_energy);
  double EAcc = (data.lepton.e() - lepton_mass) *
                    ((incident_energy - lepton_mass - MA) / (data.E - lepton_mass - MA))
                + lepton_mass;
  double Pt = data.lepton.perp();
  double P = sqrt(EAcc * EAcc - lepton_mass * lepton_mass);
  if (method_ == DarkBremMethod::ForwardOnly) {
    unsigned int i = 0;
    while (Pt * Pt + lepton_mass * lepton_mass > EAcc * EAcc) {
      // Skip events until the transverse energy is less than the total energy.
      i++;
      data = sample(incident_energy);
      EAcc = (data.lepton.e() - lepton_mass) *
                 ((incident_energy - lepton_mass - MA) / (data.E - lepton_mass - MA))
             + lepton_mass;
      Pt = data.lepton.perp();
      P = sqrt(EAcc * EAcc - lepton_mass * lepton_mass);

      if (i > maxIterations_) {
        std::cerr
            << "Could not produce a realistic vertex with library energy "
            << data.lepton.e() << " GeV.\n"
            << "Consider expanding your libary of A' vertices to include a "
               "beam energy closer to "
            << incident_energy << " GeV."
            << std::endl;
        break;
      }
    }
  } else if (method_ == DarkBremMethod::CMScaling) {
    CLHEP::HepLorentzVector el(data.lepton.px(), data.lepton.py(), data.lepton.pz(),
                               data.lepton.e());
    double ediff = data.E - incident_energy;
    CLHEP::HepLorentzVector newcm(data.centerMomentum.px(), data.centerMomentum.py(),
                                  data.centerMomentum.pz() - ediff,
                                  data.centerMomentum.e() - ediff);
    el.boost(-1. * data.centerMomentum.boostVector());
    el.boost(newcm.boostVector());
    double newE = (data.lepton.e() - lepton_mass) *
                      ((incident_energy - lepton_mass - MA) / (data.E - lepton_mass - MA)) +
                  lepton_mass;
    el.setE(newE);
    EAcc = el.e();
    Pt = el.perp();
    P = el.vect().mag();
  } else if (method_ == DarkBremMethod::Undefined) {
    EAcc = data.lepton.e();
    P = sqrt(EAcc * EAcc - lepton_mass * lepton_mass);
    Pt = data.lepton.perp();
  }

  // outgoing lepton momentum
  G4double PhiAcc = G4UniformRand()*2*pi;
  G4double recoilMag = sqrt(EAcc * EAcc - lepton_mass*lepton_mass)*GeV;
  G4ThreeVector recoil;
  double ThetaAcc = std::asin(Pt / P);
  recoil.set(std::sin(ThetaAcc) * std::cos(PhiAcc),
                             std::sin(ThetaAcc) * std::sin(PhiAcc),
                             std::cos(ThetaAcc));
  recoil.setMag(recoilMag);
  return recoil;
}

void G4DarkBreMModel::GenerateChange(
    G4ParticleChange &particleChange, const G4Track &track,
    const G4Step &step) {
  // mass of incident lepton 
  double Ml = track.GetDefinition()->GetPDGMass() / CLHEP::GeV;

  // convert to energy units in LHE files [GeV]
  G4double incidentEnergy = step.GetPostStepPoint()->GetTotalEnergy()/CLHEP::GeV;

  G4ThreeVector recoilMomentum = scale(incidentEnergy, Ml);
  recoilMomentum.rotateUz(track.GetMomentumDirection());

  // create g4dynamicparticle object for the dark photon.
  // define its 3-momentum so we conserve 3-momentum with primary and recoil
  // lepton NOTE: does _not_ take nucleus recoil into account
  G4ThreeVector darkPhotonMomentum =
      track.GetMomentum() - recoilMomentum;
  G4DynamicParticle *dphoton =
      new G4DynamicParticle(G4APrime::APrime(), darkPhotonMomentum);

  // stop tracking and create new secondary instead of primary
  if (alwaysCreateNewLepton_) {
    /*
     * Create a new lepton in order to make it easier to extract 
     * the outgoing sim-level dark brem kinematics.
     */
    G4DynamicParticle *el = new G4DynamicParticle(
        track.GetDefinition(), recoilMomentum);
    particleChange.SetNumberOfSecondaries(2);
    particleChange.AddSecondary(dphoton);
    particleChange.AddSecondary(el);
    particleChange.ProposeTrackStatus(fStopAndKill);
    // continue tracking
  } else {
    /*
     * just have primary lose energy (don't rename to different track ID)
     *
     * This branch is untested and so we are not sure if it works as expected.
     */
    particleChange.SetNumberOfSecondaries(1);
    particleChange.AddSecondary(dphoton);
    particleChange.ProposeMomentumDirection(recoilMomentum.unit());
    double recoil_energy = sqrt(recoilMomentum.mag2()+Ml*Ml);
    // energy of primary recoiling
    G4double finalKE = recoil_energy - Ml;
    particleChange.ProposeEnergy(finalKE);
  }
}

void G4DarkBreMModel::SetMadGraphDataLibrary(const std::string& path) {
  /*
   * print status to user so they know what's happening
   */
  if (GetVerboseLevel() > 0) G4cout << "[ G4DarkBreMModel ] : loading event librariy..." << G4endl;

  parseLibrary(path, aprime_lhe_id_, madGraphData_);

  if (madGraphData_.size() == 0) {
    throw std::runtime_error("BadConf : Unable to find any library entries at '"+path+"'\n"
        "  The library is either a single CSV file or a directory of LHE files.\n"
        "  Any individual file can be compressed with `gzip`.\n"
        "  This means the valid extensions are '.lhe', '.lhe.gz', '.csv', and '.csv.gz'");
  }

  MakePlaceholders();  // Setup the placeholder offsets for getting data.

  if (GetVerboseLevel() > 0) G4cout << "[ G4DarkBreMModel ] : done" << G4endl;

  /*
   * Print out loaded MG library
   */
  if (GetVerboseLevel() > 1) {
    G4cout << "MadGraph Library of Dark Brem Events:\n";
    for (const auto &kV : madGraphData_) {
      G4cout << "\t" << kV.first << " GeV Beam -> "
                << kV.second.size() << " Events\n";
    }
    G4cout << G4endl;
  }

  return;
}

void G4DarkBreMModel::MakePlaceholders() {
  currentDataPoints_.clear();
  maxIterations_ = 10000;
  for (const auto &iter : madGraphData_) {
    currentDataPoints_[iter.first] = int(G4UniformRand() * iter.second.size());
    if (iter.second.size() < maxIterations_)
      maxIterations_ = iter.second.size();
  }
}

OutgoingKinematics
G4DarkBreMModel::sample(double incident_energy) {
  // Cycle through imported beam energies until the closest one above is found,
  // or the max is reached.
  double samplingE = 0.;
  for (const auto &keyVal : currentDataPoints_) {
    samplingE = keyVal.first;  // move samplingE up
    // check if went under the sampling energy
    //  the map is sorted by key, so we can be done right after E0 goes under
    //  samplingE
    if (incident_energy < samplingE) break;
  }
  // now samplingE is the closest energy above E0 or the maximum energy imported
  // from mad graph

  // Need to loop around if we hit the end, in case our random
  // starting position happens to be late enough in the file
  if (currentDataPoints_.at(samplingE) >= madGraphData_.at(samplingE).size()) {
    currentDataPoints_[samplingE] = 0;
  }

  // increment the current index _after_ getting its entry from
  // the in-memory library
  return madGraphData_.at(samplingE).at(currentDataPoints_[samplingE]++);
}

}  // namespace g4db

