# GRASPy
A python wrapper for atomic structure calculations with GRASP 2018.

This package provides a Python interface to the GRASP 2018 program, which implements multiconfiguration Dirac-Hartree-Fock + Configuration Interaction calculations on single-atom systems.
The goal is to make complex GRASP calculations easier for the newcomer to implement.
## Why Use GRASPy?
Relativistic atomic structure calculations are both complex and computationally intensive. At many points in the multi-step calculationThere are many ways to reduce the many-body problem of interest to a tractable computation.

We hope this package will facilitate automatic exploration of the computational parameter space using tools in python like `itertools` in order to reach acceptable consistency with experimental spectroscopy data, which is especially difficult for heavy, many-electron atomic systems and requires much trial and error.

We hope that condensing complex GRASP calculations into single Python scripts will make it easy to share calculation results with others as part of publications, and make these difficult, complex calculations fully reproducible and transparent to other members of the community.

## How To Install
With a working GRASP 2018 installation and Python 3.6+, GRASPy should work right out of the box. Interfacing from GRASPy to GRASP is done entirely via the `os.subprocess()`, and `pandas.read_csv` is used to read out GRASP output (e.g. subshell energies) for manipulation and plotting in Python scripts.

## Architecture
Instead of executing a bunch of shell scripts in a particular folder as is documented in the GRASP 2018 manual,
we introduce an object-oriented framework for inputting parameters to GRASP routines.

Each routine offered in GRASP 2018 is an object, which takes in input files, some user-specified parameters, and returns some output files.
Creation of an object allows users to specify calculation parameters which were previously inputted through a config file for each script execution.
At runtime, GRASPy checks that all input files are present, and then runs the calculation in a specified working directory, printing output to the terminal--just as if you were running GRASP from the command line manually.

## Current State
So far, the following GRASP routines have working interfaces for default optionsets, and some work is being done to implement non-default calculation settings for many of these routines.

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
```
In addition, an `initialize()` function is provided to create new working directories. A sample calculation has been provided in `singlet-triplet.py`, in which we calculate the 3P1 and 1P1 levels of neutral ytterbium, demonstrating how a complex calculation with many different orbitals, each with their own DHF procedure, and a final CI run can be consolidated into a single `.py` script.

## Acknowledgments
This program was developed by Calvin Leung, with significant input from Alex Ozdemir and Joonseok Hur. If you use this tool in your work, please mention it. Feel free to contact us with any questions or requests.

