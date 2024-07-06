# -*- coding: utf-8 -*-
"""
Created on Thur Jul 04 17:27:53 2024

@author: eliot
"""


from scipy.signal import freqs
import numpy as np
import sympy as sp
from symbolicspice import Circuit, plotTransfertFunction

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

# example with the "tone" stage of the MXR Distortion +
inputList  =   ['Vin 1 0 1',
                'C1 1 0 1e-9',
                'C2 1 2 10e-9',
                'R1 2 3 10e3',
                'R2 3 0 1e6',
                'OP1 3 4 5',
                'R3 4 5 1e6',
                'C3 4 5 10e-12',
                'C4 4 6 47e-9',
                'R4 6 7 4.7e3',
                'R5 7 0 1']

#declare a circuit object
circuit = Circuit(inputList)

circuit.display_components()
circuit.stamp_system()
circuit.print_system()

'''
#Solve the system
circuit.solve_system(simplify=False, use_symengine=False)

print('\n\nx solution vector:\n')
sp.pprint(sp.Eq(circuit.x, circuit.x_solution))
'''

# Get the symbolic transfer function between the output node and the input node
# DONT FORGET TO CHANGE THE OUTPUT AND INPUT NODES ACCORDING TO YOUR NETLIST !
H = circuit.get_symbolic_transfert_function(output_node = 5, input_node = 1)

'''
print('\n\nSymbolic transfer function:\n')
sp.pprint(H.sympyExpr)
'''

b, a = H.symbolic_analog_filter_coefficients()

print(f'\nb coefficients :{b}')
print(f'a coefficients :{a}')


# Get numerical coefficients for the range of R4 values
component_values = {'C3': np.array([10e-12, 22e-12, 50e-12]),
                    'R4': np.array([4.7e3, 10e3, 22e3, 47e3, 100e3, 500e3])}
b_num, a_num =  H.numerical_analog_filter_coefficients(component_values)

f = np.arange(1, 20e3)
w = 2*np.pi*f

h = []
for i in range(len(b_num)):
    h.append([freqs(b_num[i][j], a_num[i][j], worN=w)[1] for j in range(len(a_num[i]))])
h = np.array(h)


plotTransfertFunction(f, h, semilogx=True, dB=True, phase=True)

