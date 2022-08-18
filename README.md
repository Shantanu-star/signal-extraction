# Alice Heavy flavor class for invariant mass fitting
This class can be used for fitting an invariant mass histogram. It has multiple options for signal and background functions.

## Installation
For the installion of the class kindly follow the following steps:
1. Clone this repository
- `git clone git@git.cbm.gsi.de:sh.khan/signal_extraction.git`
2. Source your ROOT
- `source /path/to/root/install/bin/thisroot.sh`
3. Make the build directory and go there 
- `mkdir build`
- `cd build`
4. Compile and install the class
- `cmake ./`
- `make`
5. 
For MAC users
uncomment line 9 in signal_extraction/src/CMakeList.txt

## Usage for the first time
To use the class, one needs to export it to root as a library

- `export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/the/path/to/build/src/`

It has a test.cpp file and one can use it to check whether everything has been done properly

- `root -l test.cpp`

## For expert user (good practices)
For a specific particle if higher order polynomial option is desired then use the following option because the high order polynomial function is defined differently than the pol2
- `SetParticlePdgMass(Particle_PDG_mass)`
