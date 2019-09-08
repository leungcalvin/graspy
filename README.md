# graspy
A python wrapper for GRASP 2018 atomic structure calculations.

This package provides a Python interface to the GRASP 2018 program. 
The goal is to make complex GRASP calculations easier for the newcomer to implement.
## Why Graspy?
Scientifically, we hope this package will facilitate automate exploration of the computational parameter space with tools in python like `itertools` in order to reach acceptable agreement with experimental data, which is especially difficult for heavy, many-electron atomic systems.

We hope that condensing complex GRASP calculations into single files will make it easy to share calculation results with others as part of e.g. publications, and make these difficult, complex calculations fully reproducible and transparent.

## Architecture
Instead of executing a bunch of shell scripts in a particular folder as is documented in the GRASP 2018 manual,
we introduce an object-oriented framework for inputting parameters to GRASP routines.

Each routine is an object, which takes in input files, some user-specified parameters, and returns some output files.
Creation of an object allows users to specify calculation parameters which were previously inputted through a config file for each script execution.
Graspy checks that all input files are present at runtime, and then runs the calculation, printing output to the terminal--just as if you were running GRASP from the command line.
When you want to execute the calculation, a working directory can be specified in order to keep intermediate results of exploratory calculations separate. 

## Current State
So far, the following GRASP routines have working interfaces for default optionsets, and some work is being done to implement non-default calculation settings for many of these routines, e.g. the zero- and first-order space formalism in `rangular`.

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
```
In addition, an `initialize()` function is provided to create new working directories. A sample calculation has been provided in `singlet-triplet.py`, demonstrating how a formerly complex calculation with many GRASP steps and many shell scripts can be consolidated into a single `.py` file.

## Acknowledgments
This program was developed by Calvin Leung, with significant input from Alex Ozdemir and Joonseok Hur. If the community likes this tool, we might end up writing a small paper and would appreciate citations. But for now, just mention `graspy` in your papers. :)

