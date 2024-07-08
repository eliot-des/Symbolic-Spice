# -*- coding: utf-8 -*-
"""
Created on Thur Jul 04 17:27:53 2024

@author: eliot
"""


from scipy.signal import freqs
import numpy as np
import sympy as sp
from symbolicspice import Circuit, plotTransfertFunction, loadnet

#declaration of the circuit
'''
# example with a LC filter of cutoff frequency 520 Hz
inputList  =   ['Vin 1 0 1',
                'L1 1 2 2.29e-3',
                'R1 2 0 8']

'''
# example with an active low-shelf filter with OPA 
inputList  =   ['Vin 1 0 1',
                'OP1 1 2 3',
                'R1 2 3 2370',
                'C1 2 3 69.5e-9',
                'R2 2 0 2359',
                'R3 3 0 8']

# example with the "tone" stage architecture of the DOD 250 overdrive pedal
# Values depends on the schematic found on the internet (but topology is correct)
'''
inputList  =   ['Vin 1 0 1',
                'C1 1 0 0.01e-9',
                'C2 1 2 10e-9',
                'R1 2 3 10e3',
                'R2 3 0 470e3',
                'OP1 3 4 5',
                'R3 4 5 1e6',
                'C3 4 5 25e-12',
                'C4 4 6 47e-9',
                'R4 6 7 4.7e3',
                'R5 7 0 1']
'''
inputList   =  ['Vin 1 0 1',
                'C1 1 2 33e-9',
                'C2 2 3 33e-9',
                'R1 3 0 4824',
                'OP1 3 4 4',
                'R2 2 4 4824',
                'R3 4 0 8',]       

#declare a circuit object
circuit = Circuit(inputList)
#circuit = Circuit(r'FMV Tone Stack.net')
#circuit.display_components()
circuit.stamp_system()
circuit.print_system()

'''
# Solve the system :
# The Code in this block can be useful if you want 
# to see the symbolic solution of the x vector.

circuit.solve_system(simplify=True, use_symengine=False)

print('\n\nx solution vector:\n')
sp.pprint(sp.Eq(circuit.x, circuit.x_solution))
'''

# Get the symbolic transfer function between the output node and the input node
# DONT FORGET TO CHANGE THE OUTPUT AND INPUT NODES ACCORDING TO YOUR NETLIST !
H = circuit.get_symbolic_transfert_function(output_node = 4, input_node = 1, normalize=True)

'''
print('\n\nSymbolic transfer function:\n')
sp.pprint(H.sympyExpr)
'''

b, a = H.symbolic_analog_filter_coefficients()

print(f'\nb coefficients :{b}')
print(f'a coefficients :{a}')


f = np.geomspace(20,20000, num=1000)
w = 2*np.pi*f


# Get numerical coefficients for the range of R4 values and C3 values.
# You can also not given a component_values argument to the function.
# In this case, the function will take the default values of the components set in the circuit object.
#3D Case :
'''
component_values = {'C3': np.array([10e-12, 22e-12, 50e-12]),
                    'R4': np.array([4.7e3, 10e3, 22e3, 47e3, 100e3, 500e3])}
b_num, a_num =  H.numerical_analog_filter_coefficients(component_values)

h = np.array([[freqs(b_num[i][j], a_num[i][j], worN=w)[1] for j in range(len(a_num[i]))] for i in range(len(b_num))])
plotTransfertFunction(f, h, legend = component_values, semilogx=True, dB=True, phase=True)
'''
'''

#2D Case
component_values = {'R1': np.array([4.7e3, 10e3, 22e3, 47e3, 100e3, 220e3]), 'R2': np.array([4.7e3, 10e3, 22e3, 47e3, 100e3, 220e3])}
b_num, a_num =  H.numerical_analog_filter_coefficients(component_values, combination='sequential')

h = np.array([freqs(b_num[i], a_num[i], worN=w)[1] for i in range(len(a_num))])
plotTransfertFunction(f, h, legend = component_values, semilogx=True, dB=True, phase=True)

'''
#1D Case
b_num, a_num =  H.numerical_analog_filter_coefficients()

_, h = freqs(b_num, a_num, worN=w)
plotTransfertFunction(f, h, legend='test', semilogx=True, dB=True, phase=True)
'''