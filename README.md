# GRASPy
A python wrapper for atomic structure calculations with GRASP 2018.

This package provides a Python interface to the GRASP 2018 program, which implements multiconfiguration Dirac-Hartree-Fock + Configuration Interaction calculations on single-atom systems.
The goal is to make long GRASP 2018 calculations easier for the newcomer and expert to implement.
## Why Use GRASPy?
Relativistic atomic structure calculations are both complex and computationally intensive. GRASP 2018's program flow using input/output files and shell scripts to program calculation runs is straightforward but error-prone as calculations get more intricate. At many points in the multi-step calculation, trying different parameter sets is required to achieve concordance with experiment. Notably, it is difficult to find a good multireference to start calculations. 


We hope that condensing complex GRASP calculations into single Python scripts will make it easy to share calculation results with others as part of publications, and make these difficult, complex calculations fully reproducible and transparent to other members of the community. In the long run, integration with packages in python like `itertools` that can automate the trial-and-error can be used reach acceptable consistency with experimental spectroscopy data, which is especially difficult for heavy, many-electron atomic systems and requires much trial and error.


## How To Install
With a working GRASP 2018 installation and Python 3.6+, GRASPy should work right out of the box. Interfacing from GRASPy to GRASP is done entirely via the `os.subprocess()`, and `pandas.read_csv` is used to read out GRASP output (e.g. subshell energies) for manipulation and plotting in Python scripts.

## How to Use
In the GRASP 2018 manual, several calculation examples are provided. As a pedagogical example, we have transcribed the first example: the calculation of the 1s^2 2s ^2 S$ and 1s^2 2p .^2P levels in Li I|into a Python script (`example1.py`). 

## Architecture
Instead of executing a bunch of shell scripts in a particular folder as is documented in the GRASP 2018 manual,
we introduce an object-oriented framework for inputting parameters to GRASP routines.

Each routine offered in GRASP 2018 is an object, which takes in input files, some user-specified parameters, and returns some output files.
Creation of an object allows users to specify calculation parameters which were previously inputted through a config file for each script execution.
At runtime, GRASPy checks that all input files are present, and then runs the calculation in a specified working directory, printing output to the terminal--just as if you were running GRASP from the command line manually.

## Current State
So far, the following GRASP routines can be called in GRASPy, and some work is being done to implement non-default calculation settings for many of these routines.

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
This program was developed by Calvin Leung, with significant input from Alex Ozdemir and Joonseok Hur. If you use this tool in your work, please mention it. Feel free to contact us with any questions or requests.

