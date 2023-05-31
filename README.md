# TRHEPD-OPT
This is a software repository for our preprint article entitled "A fast and accurate computation method for reflective diffraction simulations."

There are four folders:
* str-org: original sim-trhepd-rheed.
* str-lapack: replace some subroutined and do-loops with LAPACK's and BLAS's one, and parallelize some part o code by OpenMP.
*  sim-trhepd-rheed: new one.
* test: all the test results and scripts.

## Requirements
You need a fortran compiler and LAPACK implementation for compilation.
A python environment is needed for running script files in the plot folder.

## LICENSE
See the LICENSE file for each folder. Scripts under test folder are public domain.
We thank to Izumi Motizuki at KEK for giving us the input files for benchmarking test listed below:
* test/test23/bulk.txt
* test/test23/surf.txt
* test/test47/bulk.txt
* test/test47/surf.txt
* test/test521/bulk.txt
* test/test521/surf.txt
