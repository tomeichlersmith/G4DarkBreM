# G4DarkBreM

## Installation
The only dependency of G4DarkBreM is Geant4 which has [an extensive installation guide](http://cern.ch/geant4-userdoc/UsersGuides/InstallationGuide/html/). While G4DarkBreM is not explicitly limited to a certain platform, it has only been used on Linux-based operating systems and explicitly uses C++11 standard features.

After installing Geant4, one can build and install G4DarkBreM using tools probably used to install Geant4 (if built from scratch).
```
cmake -B build -S . -DCMAKE_INSTALL_PREFIX=<my-install>
cd build
make install
```

Additionally, G4DarkBreM can be pulled into a more expansive simulation framework as a submodule and included as a subdirectory in CMake 
```
add_subdirectory(G4DarkBreM)
```
This defines the `G4DarkBreM` cmake target which later targets can link to, for example
```
target_link_libraries(MySim PUBLIC G4DarkBreM)
```

## General Comments
This dark brem simulation method functions by writing a `G4VDiscreteProcess`.
The brem process `G4DarkBremmstrahlung` has two important handles into the simulation
that it inherits from `G4VDiscreteProcess`.

1. `GetMeanFreePath` --- calculate an effective mean free path given the current state of the particle
2. `PostStepDoIt` --- actually perform the dark brem interaction

These are the two core additional methods; and the ones that require the most complexity. 
In order to separate out these core, complex methods from the other process methods that are required to interface with Geant4, 
the cross section calculations and dark brem simulation is done in a "Dark Brem Model" which is owned by the process. 
The class `G4DarkBremsstrahlungModel` is an abstract class detailing the necessary methods for defining a specific dark brem model. 
One specific model is defined.

1. `G4DarkBreMModel` -- the new model using an imported library of dark brem events as described in the paper

This design choice was implemented with the assumption that future, alternative dark brem models can be implemented.
Additionally, a particle class `G4APrime` is defined to allow the use of the A' particle in Geant4.
Additional downstream classes (biasing and filtering) are implemented in other modules since they are detector-specific.

In the table below, the configuration options are given. Some options (especially ones specific to the model),
require more elaboration which is given in sections below.

Parameter           | Required? | Level   | Description
--------------------|-----------|---------|------------
MadGraph Data       | Yes       | Model   | Path to directory containing LHE files of dark brem events
A' Mass             | Yes       | Process | Mass in MeV of the A' to simulate (must match mass in LHE data)
Scaling Method      | No        | Model   | How the model should interpret the LHE data (more detail below; Default = Forward Only)
Mixing Strength     | No        | Model   | Cross section mixing strength epsilon (Default = 1)
Threshold           | No        | Model   | Minimum incident energy in GeV for cross section to be non-zero (Default = twice A' mass)
Allow >1 Dark Brem  | No        | Process | Allow for dark brem process to occur more than once (Default = False)
Cache Cross Section | No        | Process | Cache the calculated cross section (Default = True)

### Cross Section Calculation
The estimate for the total cross section given the material and the lepton's energy is done using the WW approximation as mentioned earlier. 
The WW approximation is done using Boost's ODE library to numerically calculate the integrals. The actual formulas are listed here for reference.

$$
\sigma = 4 \frac{pb}{GeV}\alpha_{EW}^3 \int^{m_A^2}_{m_A^4/(4E_0^2)} \chi(t)dt \int_0^{\min(1-m_e/E_0,1-m_A/E_0)} \frac{d\sigma}{dx}(x)dx
$$

where

$$
\chi(t) = \left( \frac{Z^2a^4t^2}{(1+a^2t)^2(1+t/d)^2}+\frac{Za_p^4t^2}{(1+a_p^2t)^2(1+t/0.71)^8}\left(\frac{1+t(m_{up}^2-1)}{4m_p^2}\right)^2\right)\frac{t-m_A^4/(4E_0^2)}{t^2}
$$

$$
a = \frac{111.0}{m_e Z^{1/3}}
\quad
a_p = \frac{773.0}{m_e Z^{2/3}}
\quad
d = \frac{0.164}{A^{2/3}}
$$

$$
\frac{d\sigma}{dx}(x) = \sqrt{1-\frac{m_A^2}{E_0^2}}\frac{1-x+x^2/3}{m_A^2(1-x)/x+m_e^2x}
$$

$E_0$ is the incoming lepton's energy in GeV, 
$m_e$ is the mass of the lepton in GeV, 
$m_A$ is the mass of the dark photon in GeV, 
$m_p = 0.938$ is the mass of the proton in GeV, 
$m_{up} = 2.79$ is the mass of the up-quark in MeV, 
$A$ is the atomic mass of the target nucleus in amu, 
$Z$ is the atomic number of the target nucleus, 
$\alpha_{EW} = 1/137$ is the fine-structure constant, 
and $pb/GeV = 3.894\times10^8$ is a conversion factor from GeV to pico-barns.

This approximation is not perfect; however, it follows the trend expected from generating events in MadGraph. 

### MadGraph Data
The reason the process allows for a directory of LHE files is so that multiple different incident energy leptons can be given samples of events. 
This allows the process to always choose a sampled energy near to the actual current energy of the lepton. 
Keeping the difference between the actual energy and the MadGraph incident energy close helps keep all of the scaling methods described below more accurate.

### Scaling Method
Three options for different interpretation or "scaling" methods are implemented.

1. Forward Only : Loops through the LHE data until the transverse momentum squared and the mass squared from the LHE total less than the actual energy squared. Uses this chosen transverse momentum with the actual energy choosing the longitudinal momentum to be positive.
2. CM Scaling : Scales the LHE vertex momenta to the actual lepton energy using Lorentz boosts into and out-of the Center of Momentum frame.
3. Undefined : Just uses the LHE data without any modification (i.e. assumes the actual lepton 4-momentum is the 4-momentum in the LHE file)

While energy fraction and transverse momentum are found to be accurate scaling variables, 
they do not define the sign of the outgoing longitudinal momentum. 
To solve this, the simulation always assumes that the lepton goes forwards and this scaling method is referred to as "Forward Only". 
To recover these backwards scattering events, an alternate scaling method has been developed using the center of mass of the lepton/dark photon system. 
Study of MadGraph events with varying incident lepton energies showed that the transverse momentum of the center of mass of the scattered lepton and 
dark photon changes very slowly, and that its longitudinal momentum scaled linearly with incident energy. 
Applying this information, the "CM Scaling" method works by finding the center of mass of the lepton and dark photon from the sampled MadGraph event, 
then reducing its longitudinal momentum and total energy by the difference between the sampled and desired incident energies. 
The sampled lepton is then boosted from the initial center of mass to the new one, and its outgoing angle in the new frame is collected. 
While using this angle preserves backwards-going events, the total energy was less reliable than using the initial method. 
To reduce this, after finding the outgoing angle of the scattered lepton the total energy is sampled by preserving the fraction of kinetic energy 
in the same way as the forward only method. Even with this modification, 
the cost of keeping the backwards-going leptons is seen in slightly less accurate longitudinal momenta.

All of these methods have their advantages and dis-advantages, so they are all left in as options. 
In all cases, using a fine-grained incident lepton energy library is **heavily suggested** in order to minimize the discrepancies arising from any of these scaling methods.
