# GRASPy
A python wrapper for atomic structure calculations with GRASP 2018.

This package provides a Python interface to the GRASP 2018 program, which implements multiconfiguration Dirac-Hartree-Fock + Configuration Interaction calculations on single-atom systems.
The goal is to make long GRASP 2018 calculations easier for the newcomer and expert to implement.
## Why Use GRASPy?
Relativistic atomic structure calculations are both complex and computationally intensive. GRASP 2018's program flow using input/output files and shell scripts to program calculation runs is straightforward but error-prone as calculations get more intricate. At many points in the multi-step calculation, trying different parameter sets is required to achieve concordance with experiment. Notably, it is difficult to find a good multireference to start calculations. 


We hope that condensing complex GRASP calculations into single Python scripts will make it easy to share calculation results with others as part of publications, and make these difficult, complex calculations fully reproducible and transparent to other members of the community. In the long run, integration with packages in python like `itertools` that can automate the trial-and-error can be used reach acceptable consistency with experimental spectroscopy data, which is especially difficult for heavy, many-electron atomic systems and requires much trial and error.


## How To Install
With a working GRASP 2018 installation and Python 3.6+, GRASPy should work right out of the box. Interfacing from GRASPy to GRASP is done entirely via the `os.subprocess()`, and `pandas.read_csv` is used to read out GRASP output (e.g. subshell energies) for manipulation and plotting in Python scripts.

```
git clone https://github.com/leungcalvin/graspy/
```

## How to Use
As a pedagogical example, we have transcribed the first example , which calculates the 1s2 2s 2S and 1s2 2p 2P levels in Li I, into a Python script (`example1.py`). Running `python example1.py` will produce a directory `example1/` and output files inside.

## Architecture
Instead of executing a bunch of shell scripts in a particular folder as is documented in the GRASP 2018 manual,
we introduce an object-oriented framework for inputting parameters to GRASP routines.

Each routine offered in GRASP 2018 is an object, which takes in input files, some user-specified parameters, and returns some output files.
Creation of an object allows users to specify calculation parameters which were previously inputted through a config file for each script execution.
At runtime, GRASPy checks that all input files are present, and then runs the calculation in a specified working directory, printing output to the terminal--just as if you were running GRASP from the command line manually.

## Current State
So far, the following GRASP routines can be called in GRASPy. Not all non-default options have been implemented.

```
rnucleus
rcsfgenerate
rcsfinteract
rangular
rwfnestimate
rmcdhf
rci
rmixextract
jj2lsj
rlevels
rcsfzerofirst
rmixaccumulate
rbiotransform
rhfs
rtransition
```

## Acknowledgments
This core of this program was developed by Calvin Leung (myself), with significant input from Alex Ozdemir and Joonseok Hur. Feel free to contact me with any questions, feedback, or special requests.

