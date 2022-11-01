# g4db executables

G4DarkBreM comes with three executables that are helpful for studying its behavior.

All of these executables output CSV text files for easier analysis.

## g4db-scale
This executable calls the sample and scale procedure _directly_, allowing the user to study how the procedure affects the outgoing kinematic distributions without having to wait for an entire Geant4 simulation to progress.

## g4db-xsec-calc
This executable, similar to above, allows the user to call the cross section calculation directly so that the user can validate and test the calculation.

## g4db-simulate
This is a full Geant4 simulation focused on a simple prism of material limited to electrons or muons shot directly into it. This is not G4DarkBreM's only use case, but it is a good one for testing that it is functioning properly.

This simulation is over simplified and **should not** be used for production-level studies. Among other simplifications, it does not implement a biasing procedure for the dark brem process and does not alter the random number seed.
