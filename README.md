# AMS02AH3

### Summary
This code, consisting of the script ah3_mc.C, was created to simulate the coalescence of antihelium-3 nuclei in the decay of W-bosons. Using information on Z-boson decay from the Particle Data Group, and coalescence momentum information from Carlson et al (arXiv:1401.2461v2), it makes a number of assumptions to provide a rough model for antihelium-3 coalescence.

### Usage
The script is meant to be run by calling the function ConvergeMulti in ROOT. Below is the syntax:
```sh
void ConvergeMulti(const UInt_t n = 7, const UInt_t range = 3, const UInt_t start = 3)
```
  - n is the maximum number of events to run (per # of antinucleons produced per event)
  - range is the range of # of antinucleons per event to run over (increasing from "start")
  - start is the initial value of # of particles per event to run over (start >= 3)

### Dependencies
[![N|Solid](https://d35c7d8c.web.cern.ch/sites/d35c7d8c.web.cern.ch/files/website-banner-allnew-croped_3.png)](https://root.cern.ch)

This script makes use of ROOT, the data analysis framework developed at CERN (see https://root.cern.ch for information on installing ROOT).
The code was run in the ROOT interpreter. For example, one may run the following:
```sh
$ root -l
$ root [0] .L ah3_mc.C
$ root [1] ConvergeMulti(6, 2, 4)
$ root [2] GenHist("50k_events_with_errors.root")
```
One may alternatively call Converge(n, m) to simply run a maximum of n events, where each event produces m antinucleons. 

### Acknowledgements

I would like to thank Professor Daniel Marlow for serving as my advisor for this project. I would also like to thank my father Edmond Offermann for teaching me most of what I know about ROOT's functionality, and Dr. Ryan R. Rios of NASA Johnson Space Center for teaching me good ROOT coding practices.

### Contact
For any questions, I may be reached at jano@princeton.edu. Thank you.

Jan Tuzlic Offermann
