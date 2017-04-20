# README #

This is a multi-objective integer programming (MOIP) algorithm, which is heavily based on the work discussed in the following paper.

Ozlen, M., Burton, B. A., & MacRae, C. A. (2014). Multi-objective integer programming: An improved recursive algorithm. Journal of Optimization Theory and Applications, 160(2), 470-482.

The original algorithm implementation can be found at https://bitbucket.org/melihozlen/moip_aira/ but this version has significant changes, noticeably a switch to C++ and boost for various utility functions, as well as the use of parallel processing.

### What is this repository for? ###

This software can be used to optimise various multi-objective integer programming problems. Currently, work is not yet stabilised and no guarantees can be offered on the correctness of various versions of the implementation. Official releases will hopefully be coming shortly.

### How do I get set up? ###

The implementation works on Linux operating systems, and requires IBM ILOG CPLEX (tested with versions 12.6.3 and 12.7.0) and Boost (specifically the program_options library). It also uses CMake (version 3.4 or higher) for building.

#### Compiling

Clone the repository

    $ git clone https://github.com/WPettersson/moip_aira
    
Create a directory in which to build

    $ mkdir build && cd build
    
Set up the build environment

    $ cmake ../ -DCPLEX_ROOT=/path/to/IBM/ILOG/CPLEX_Studio127/
    
Compile

    $ make

Your executable now resides in the `src` directory. If you wish, you can run `make install` to install it to `/usr/local/bin`, assuming you have rights to do so. You can also run the included test suite by calling `make test`.

moip_aira uses an extended LP file format where multiple objectives are defined as additional constraints after the original problem's constraints. The right-hand-side value of the last constraint defines the number of objectives. Example LP files are provided under a  separate folder.

### Who do I talk to? ###

Dr William Pettersson (william@ewpettersson.se) is the lead developer of this particular implementation of this algorithm. For more details on the original algorithm, you may also wish to talk to Assoc. Prof. Melih Ozlen (melih.ozlen@rmit.edu.au)
