
ChemPotPy, CHEMical library of POTential energy surfaces in PYthon 
==================================================================

Aug. 26, 2023

Authors: Yinan Shu, Zoltan Varga, Dayou Zhang, Donald G. Truhlar
University of Minnesota, Minnesota, United States

ChemPotPy is a library for analytic representation of single-state 
and multi-state potential energy surfaces and couplings. 

All fortran source code are stored in folder chempotpy 


How to install
--------------

* The users can either compile all .so modules by yourself:

  1. all potential energy surface subroutines are located in chempotpy/chempotpy/
  2. to compile each of these surface subroutines yourself, execute install.script

* Or, the users can use pre-compiled all .so modules, Follow the following steps:

  1. install Conda 
  2. create a virtural environment that uses latest python and gfortran,
     and call it, for example, chempotpy
     "conda create -n chempotpy"
  3. install gfortran 
     "conda install -c conda-forge gfortran"
  4. install numpy 
     "conda install numpy"
  5. "pip install chempotpy"


Citation
--------

The following paper should be cited in publications utilizing the
ChemPotPy library in addition to the original paper that publishes 
the potential energy surface subroutine:

Shu, Y.; Varga, Z.; Truhlar, D. G.
"ChemPotPy: A Python Library for Analytic Representation of Potential 
Energy Surfaces and Diabatic Potential Energy Matrices"
to be submitted



TO CONTRIBUTE YOUR POTENTIAL
----------------------------
* Option1
send email to one of the maintainers:
Yinan Shu, yinan.shu.0728@gmail.com
Zoltan Varga, zoltan78varga@gmail.com
Dayou Zhang, zhan6350@umn.edu
Donald G. Truhlar, truhlar@umn.edu
 
* Option 2
Submit tickets on the [gitissues](https://github.com/shuyinan/chempotpy/issues)
