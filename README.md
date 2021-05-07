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

# Running a test job

The TABOO configuration is specified in the `task_*.dat` files, which correspond to the three TABOO _tasks_: a spectral analysis (`task_1.dat`), the simulation of a simple load (`task_2.dat`) and the simulation of a complex load (`task_3.dat`). A detailed description of the configuration file for each task is given in the TABOO User guide. The three default configuration files provided with TABOO set up a spectral analysis for the Earth model M3-L70-V01, which is part of the _benchmark test suite_ discussed in Spada et al. (2011). To run this analysis, it is sufficient to invoke TABOO by typing 

    $ ./taboo.exe
    
TABOO will create various output files, whose meaning is illustrated in the User guide. To plot the relaxation spectrum of the model, contained in the `spectrum.dat` output file, type the following commands:

    $ gnuplot
    gnuplot> set logscale x
    gnuplot> set logscale y
    gnuplot> plot 'spectrum.dat' using 1:5;
    
# References

Spada, G., Barletta, V., Klemann, V., Riva, R.E.M., Martinec, Z., Gasperini, P., Lund, B., Wolf, D., Vermeersen, L.L.A. and King, M., A benchmark study for glacial isostatic adjustment codes, Geophysical Journal International, 185 (1), 106-132, 2011.
