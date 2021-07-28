# README #

This is a multi-objective integer programming (MOIP) algorithm, which is heavily based on the work discussed in
Ozlen, M., Burton, B. A., & MacRae, C. A. (2014). Multi-objective integer programming: An improved recursive algorithm. Journal of Optimization Theory and Applications. doi: 10.1007/s10957-013-0364-y

The parallel elements are documented in
W. Pettersson, M. Ozlen (2019). Multi-objective integer programming: Synergistic parallel approaches. INFORMS Journal on Computing 2019; doi: 10.1287/ijoc.2018.0875 

The original algorithm implementation can be found at https://bitbucket.org/melihozlen/moip_aira/ but this version has significant changes, noticeably a switch to C++ and boost for various utility functions, as well as the use of parallel processing.

### What is this repository for? ###

This software can be used to optimise various multi-objective integer programming problems. Currently, work is not yet stabilised and no guarantees can be offered on the correctness of various versions of the implementation. Official releases will hopefully be coming shortly.

### How do I use it?

The program itself is called `aira` and is a command-line program. A complete listing of usage options can be found by simply running `aira` from a command prompt. Only one parameter is required, the program file to read. All other options tweak the parallelisation options and algorithm selection.

#### Algorithm selection

moip_aira has two in-built parallelisation algorithms, Efficient Projection Parallelisation (EPP) from https://arxiv.org/abs/1704.08417 and the synergistic parallelisation from https://arxiv.org/abs/1705.03112. The default is the synergistic parallelisation, as it performs better. Both algorithms have a few common options:

* `-t` controls how many workers the parallel algorithm should start
* `-c` controls how many threads each worker will use when calling CPLEX

Note that the above means that the total number of threads to use is `t * c`. For instance, `-t 4 -c 3` will utilise up to 12 threads at once. Both of these options default to 1.

##### EPP

EPP splits the solution space up into `t` divisions, based on the possible values the last objective function can reach. One worker is then started for each such division.

##### Synergistic

For a proper understanding of this approach, I highly recommend reading the paper https://arxiv.org/abs/1704.08417 as it is not simple to describe. This algorithm has one more possible option, `--spread`. This option selects between the "clustering" and "spreading" algorithms as described in the paper. The default is "spreading" as this was the faster algorithm in the paper, but to select "clustering" simply pass `--spread=0` to `aira`.

##### How do I specify a problem.

moip_aira uses an extended LP file format where multiple objectives are defined as additional constraints after the original problem's constraints. The right-hand-side value of the last constraint defines the number of objectives. Example LP files are provided under the examples folder.

##### No, really, just give me an example

Fine. Running `aira -p Examples/3AP05.lp` will solve the problem specified in the file `Examples/3AP05.lp`. If you inspect this file, you will note that the last constraint is of the form "... ≤ 3". This 3 indicates that the last 3 constraints in the file are not actually constraints to our problem, but are instead 3 distinct objective functions.

The output will be written to `Examples/3AP05.out`. This will show some timing and algorithm information, but mostly it will show rows of output containing 3 numbers. Each of this is a Pareto solution. Note that you can write the output to any file you choose by using the `-o MyOutput.txt` option.

### How do I install it?

The implementation works on Linux operating systems, and requires IBM ILOG CPLEX (tested with versions 12.6.3 and 12.7.0) and Boost (specifically the program_options library). It also uses CMake (version 3.4 or higher) for building.
This implementation should work just fine on Windows and OSX as well, however I do not have the ability to test those builds. If you can give me access to machines and software with which to test these, I will gladly try my best to provide some help.

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


### Who do I talk to? ###

Dr William Pettersson (william@ewpettersson.se) is the lead developer of this particular implementation of this algorithm. For more details on the original algorithm, you may also wish to talk to Assoc. Prof. Melih Ozlen (melih.ozlen@rmit.edu.au)
