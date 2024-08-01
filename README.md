# Symbolic-Spice: A python script for symbolic analysis of electronic circuits
---

The goal of Symbolic-spice is to provide a simple solution to analyse linear circuit in a symbolic fashion, in order, for example, to **derive more easily the analytical expression of a transfert function on a given circuit**. Also, Symbolic-spice can substitute (thanks to sympy) the symbolic component values by there associated numerical value, in order to plot a transfert function between two nodes/branches of this circuit.

## Under the hood of Symbolic-Spice

Symbolic-spice relies on the Modified Nodal Analysis, introduced by [Chen-When Ho and al.](https://cseweb.ucsd.edu/classes/fa04/cse245/Reading/MNA.pdf), in order to automaticaly right down the equations governing the linear electronic circuit, in a symbolic way thanks to the [sympy](https://github.com/sympy/sympy) library.

## How to use Symbolic-Spice with LT Spice
1. Create a linear schematic in LT-Spice.
2. Set a voltage source named 'Vin'.
    a. The positive side of the source needs to be connected in the desired input node.
    b. The AC stimulus of the source needs to be set to 1.
3. (Optional) Set a net label titled 'out' for setting the output as the last node.

## To do:
- [x] Add parser for LT Spice .net files.
	- [x] Add support fot V1 as default input.
    - [ ] Add support for parameters.
    - [ ] Test parser with all the components.
		- [x] Voltage source.
		- [ ] Current source.
		- [x] Resistor.
		- [x] Capacitor.
		- [x] Inductor.
		- [ ] Op-Amp.
		- [ ] Transformer.
		- [ ] Girator.
- [x] Convert scripts into library.
- [x] Add Z-Transform for analog coefficients.
    - [x] Merge into transfer function class
    - [x] Clean-up implementation
    - [x] Add normalizing.
    - [x] Solve numpy object/float situation.
- [ ] Add feature to get State-Space model from netlist.
- [ ] Add optimizer to fine tune parameters based on a target measurement.
- [ ] Use Pultec's EQP-1A as an example.
