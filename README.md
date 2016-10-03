# README #

This is a C implementation of An Improved Recursive Algorithm (AIRA) for Multi-Objective Integer Programming (MOIP) as discussed in the following paper.

Ozlen, M., Burton, B. A., & MacRae, C. A. (2014). Multi-objective integer programming: An improved recursive algorithm. Journal of Optimization Theory and Applications, 160(2), 470-482.

### What is this repository for? ###

Version: 0.9

### How do I get set up? ###

The implementation works on Linux operating systems, and requires IBM ILOG CPLEX 12.6.3 and glib-2.0.

It uses an extended LP file format where multiple objectives are defined as additional constraints after the original problem's constraints. The right-hand-side value of the last constraint defines the number of objectives. Example LP files are provided under a  separate folder.

### Who do I talk to? ###

Assoc. Prof. Melih Ozlen (melih.ozlen@rmit.edu.au)