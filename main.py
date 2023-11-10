# -*- coding: utf-8 -*-
"""
Created on Sun Oct 22 14:25:17 2023

@author: eliot
"""

from lib import create_component, display_components, create_A_matrix, create_x_vector, create_z_vector, solve_x_vector
import sympy as sp


sp.init_printing()


netlist = ['0 1 Vg', 
           '1 2 R1',
           '2 3 L1',
           '3 0 C1']


solution = solve_x_vector(netlist)

sp.pprint(solution)



#creation of the objects of the different subclasses from the netlist 
components_list = [create_component(item) for item in netlist]


#display the interpretation of the netlist by the algorithm
display_components(components_list)


#system of the form Ax = z
A = create_A_matrix(components_list)
x = create_x_vector(components_list)
z = create_z_vector(components_list)


print('\n\nA matrix :\n ')
sp.pprint(A)
print('\n\nx vector :\n ')
sp.pprint(x)
print('\n\nz vector :\n ')
sp.pprint(z)



solution =  sp.simplify(A.LUsolve(z))




sp.pprint(sp.Eq(x , solution))

#print Transfert function across R1
sp.pprint(sp.simplify((solution[0]-solution[1])/solution[0]))

