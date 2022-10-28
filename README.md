# G4DarkBreM

Geant4 Dark Bremmstrahlung from MadGraph

## Installation
The only dependency of G4DarkBreM is Geant4 which has [an extensive installation guide](http://cern.ch/geant4-userdoc/UsersGuides/InstallationGuide/html/). 
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
