# Symbolic-Spice
This python program provide a solution to obtain a symbolic analysis of your linear system thanks to the Modified Nodal Analysis.
All you have to do is to provide a netlist to the python `.main` file, and the program will return the symbolic values of the nodal voltages of the circuit.
In the netlist, you must say that the ground correspond to the node 0.

The algorithm is based on the one describe in the book frow Lawrence Nagel, called "SPICE2 : A Computer Program To Simulate Simeconductor Circuits".
If you don't want to read this book, Eric Cheever from the Swarthmore College have written a short and understandable article on this page "https://lpsa.swarthmore.edu/Systems/Electrical/mna/MNA1.html", to now how to implement this algorithm.




NOTE: this is a work in progress project.

