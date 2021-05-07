# TABOO
TABOO - a posT glAcial reBOund calculatOr

TABOO is a Fortran program for solving various problems in the context of postglacial rebound modeling. Its main functions are: i) computing loading and tidal Love numbers in normal-mode multi-exponential form; ii) simulating the time-dependent response of the Earth to a simple surface load; iii) obtaining the cumulated response of the Earth to a complex aggregate of surface loads.

# Building

TABOO does not depend on any external software package. It has been tested with recent versions of the GNU gfortran compiler, but in principle it should be possible to compile it with any Fortran compiler supporting the quadruple-precision arithmetic (`REAL*16`). While not strictly required, the [Generic Mapping Tools](https://www.generic-mapping-tools.org/) and/or [Gnuplot](http://www.gnuplot.info/) may be helpful to visualize results obtained with TABOO.

To build TABOO with the gfortran compiler, just type

    $ gfortran -O taboo.f90 -o taboo.exe
   
On Windows systems, TABOO can be compiled with gfortran either within the [Cygwin environment](http://cygwin.com) or the Windows Subsystem for Linux (WSL).

# Documentation

A user manual and a theory booklet describing the mathematical foundation of TABOO are both available inside the `doc/` directory.
