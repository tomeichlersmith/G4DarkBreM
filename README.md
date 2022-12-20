# G4DarkBreM

Geant4 Dark Bremmstrahlung from MadGraph

<p align="center">
  <a href="https://www.apache.org/licenses/LICENSE-2.0" alt="Apache 2.0 license">
    <img src="https://img.shields.io/badge/license-Apache%202-blue" />
  </a>
  <a href="https://github.com/tomeichlersmith/G4DarkBreM/actions" alt="Actions">
    <img src="https://img.shields.io/github/workflow/status/tomeichlersmith/G4DarkBreM/test.yml?branch=main" />
  </a>
  <a href="https://github.com/tomeichlersmith/G4DarkBreM/releases" alt="Releases">
    <img src="https://img.shields.io/github/v/release/tomeichlersmith/G4DarkBreM" />
  </a>
</p>

The version submitted in 2022 to Computer Physics Communications is [release v1.1.1](https://github.com/tomeichlersmith/G4DarkBreM/releases/tag/v1.1.1).
Ongoing development of this package is maintained at on GitHub [tomeichlersmith/G4DarkBreM](https://github.com/tomeichlersmith/G4DarkBreM).

## Installation
The only dependencies of G4DarkBreM are Geant4 which has [an extensive installation guide](http://cern.ch/geant4-userdoc/UsersGuides/InstallationGuide/html/)
and Boost which can be installed [from the website](https://www.boost.org/doc/libs/1_80_0/more/getting_started/unix-variants.html)
or via your package manager (e.g. [on Ubuntu](https://stackoverflow.com/questions/12578499/how-to-install-boost-on-ubuntu)).

As defined by [our CMake infrastructure](CMakeLists.txt), the minimum versions of these dependencies are 1.68 for Boost and 10.2.3 for Geant4.
We use the Boost.Math and Boost.Iostreams subcomponents of Boost if you wish to limit the size of the Boost needed to be installed.

While G4DarkBreM is not explicitly limited to a certain platform, it has only been used on Linux-based operating systems and explicitly uses C++11 standard features.

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

## Validation
Analysis and validation of G4DarkBreM has been studied in another repository 
[tomeichlersmith/ldmx-sim-technique](https://github.com/tomeichlersmith/ldmx-sim-technique).
This repository also stands as an example for integrating G4DarkBreM into a larger simulation and processing framework.
