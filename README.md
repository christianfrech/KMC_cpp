## KMC_cpp

This is a C++ and python library for both serial and parallel implementations of the BKL (rejection-free) Kinetic Monte Carlo algorithm. The lattice is specifically constructed for a BCC crystalline system, and the moves/actions in the system include vacancy-mediated diffusion along the {111} and {100} axes, as well as stripping at a pre-specified interface.

### Format of the input file 

For initializing a simulation, one must follow the input format provided in the example files. The rate catalog for all moves is determined in the header of the input file and follows the following structure:  
  
lattice_dims: [xdimension] [ydimension] [zdimension]   
atomtypes: [atom idx 0]:[atomtype name 0] [atom idx 1]:[atomtype name 10]  
num_regions: [number of regions]  
regions begin  
[region index]: [region type] [bias direction] xmin:[xmin value] xmax:[xmax value] ymin:[ymin value] ymax:[ymax value] zmin:[zmin value] zmax:[zmax value] rate_neg [bias rate along neg. direc.] rate_pos [bias rate along pos. direc.] INTERFACE(optional) 1 0 1 2 7e8  
regions end  
rates begin  
diag [{111} direction move rate] lateral [{100} direction move rate] void_threshold [number of vacs. for void diss. rate] void_rate [{void diss. rate]  
terrace_rate_111 [terrace move rate along {111} direc.] terrace_rate_100 [terrace move rate along {100} direc.] void_gb_diss_rate [void diss. rate into gb core]  
rates end  

### Organization of the project

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

