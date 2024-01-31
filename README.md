## KMC_cpp

This is a C++ and python library for a serieal-implementation of the BKL (rejection-free) Kinetic Monte Carlo algorithm. The lattice is specifically constructed for a BCC crystalline system, and the moves/actions in the system include vacancy-mediated diffusion along the {111} and {100} axes, as well as stripping at a pre-specified interface.

For initializing a simulation, one must follow the input format provided in the example files/

### Organization of the  project

The project has the following structure:

    shablona/
      |- README.md
      |- kmc_cpp/
         |- kmc_stripping.hpp/
         |- kmc_stripping_call.cpp/
	 |- hpp_files/
            |- 
	    |- 
	    |- 
         |- tests/
            |- ...
      |- doc/
         |- Makefile
         |- conf.py
         |- sphinxext/
            |- ...
         |- _static/
            |- ...
      |- setup.py
      |- .travis.yml
      |- .mailmap
      |- appveyor.yml
      |- LICENSE
      |- Makefile
      |- ipynb/
         |- ...

