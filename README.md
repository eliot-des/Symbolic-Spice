# Symbolic-Spice: A python script for symbolic analysis of electronic circuits
---

The goal of Symbolic-spice is to provide a simple solution to analyse linear circuit in a symbolic fashion, in order, for example, to derive more easily the analytical expression of a transfert function on a given circuit. Also, Symbolic-spice can substitute (thanks to sympy) the symbolic component values by there associated numerical value, in order to plot a transfert function between two nodes/branches of this circuit.

## Under the hood of Symbolic-Spice

Symbolic-spice relies on the Modified Nodal Analysis, introduced by [Chen-When Ho and al.](https://cseweb.ucsd.edu/classes/fa04/cse245/Reading/MNA.pdf), in order to automaticaly right down the equations governing the linear electronic circuit, in a symbolic way thanks to the [sympy](https://github.com/sympy/sympy) library.
