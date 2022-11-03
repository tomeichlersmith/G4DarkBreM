# Dark Brem Event Libraries

These **example** dark brem event libraries were generated using v4.3 of 
[tomeichlersmith/dark-brem-lib-gen](https://github.com/tomeichlersmith/dark-brem-lib-gen)
and then extracted to a CSV and compressed (process below).

They are merely here to help test run the simulation and be an example.
They should not be used for any production-level physics studies.

## Generation Process
Follow [Quick Start](https://github.com/tomeichlersmith/dark-brem-lib-gen#quick-start)
to setup dark brem library generation environment making sure to use `v4.2`.

Generate the electron library.
```
dbgen run --nevents 50000 \
  --max_energy 4. \
  --min_energy 0.2 \
  --apmass 0.1 \
  --target tungsten \
  --lepton electron
```

Generate the muon library.
```
dbgen run --nevents 50000 \
  --max_energy 100 \
  --min_energy 2 \
  --apmass 1 \
  --target copper \
  --lepton muon
```

The directories created by these runs can already be provided as a dark brem event library;
however, if you wish to shrink the size of the files you are carrying around, you can use
`g4db-extract-library` and `gzip` to save disk space.

First, pull out the library into a CSV file.
```
g4db-extract-library <path-to-db-lib>
```
This already saves a factor of ~10 in disk space since the scaling method doesn't use all
of the kinematic information output by MadGraph/MadEvent. Without explicitly providing
an output file, it will simply call the CSV file the same name as the directory of the 
LHE files with the `.csv` extension added. This CSV file can also be loaded by provided
as a dark brem event library.

If you wish to save another factor of ~4 in disk space, you can also compress the CSV.
```
gzip <path-to-db-lib-csv>
```
This `gzip`-compressed CSV file can also be provided as a dark brem event library.
