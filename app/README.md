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
Nevertheless, it does hold several helpful example classes in the g4db::example namespace that can be a good starting point for developing a more thorough simulation.

To view the results of an example simulation that match similar criteria as the simulation used in the paper,
you can run this executable with the provided dark brem event libraries stored in the [data](../data/) subdirectory.

For electrons
```
g4db-simulate \
  --ap-mass 0.1 \
  --depth 18 \
  --target G4_W \
  --output electrons_thick_target.csv \
  --beam 4 \
  data/electron_tungsten_MaxE_4.0_MinE_0.2_RelEStep_0.1_UndecayedAP_mA_0.1_run_3000.csv.gz \
  50000 
```

and similarly for muons
```
g4db-simulate \
  --muons \
  --ap-mass 1.0 \
  --depth 2000 \
  --target G4_Cu \
  --output muons_thick_target.csv \
  --beam 100 \
  data/muon_copper_MaxE_100.0_MinE_2.0_RelEStep_0.1_UndecayedAP_mA_1.0_run_3000.csv.gz \
  50000 
```
**Note**: For muons, the simulation will take significantly longer. This is because the cross section calculation 
for muons is significantly more complex and requires more time to calculate.

## g4db-extract-library
This helps test the library parsing procedure by reading in LHE (or `gzip` compressed LHE) into memory and then dumping the resulting library to a CSV text file. 
The output CSV can then be used by the model if the user so wishes and/or used for easier analysis of the raw library kinematics.
See g4db::parse::csv for an explanation of the columns of the CSV.
