# -*- coding: utf-8 -*-
"""
Created on Thur Jul 04 17:27:53 2024

@author: eliot
"""


from scipy.signal import freqs, freqz
import numpy as np
import sympy as sp
from symbolicspice import Circuit, plotTransfertFunction
import matplotlib.pyplot as plt

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
'''
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
H = circuit.tf(output_node = 5, input_node=1, norm=True)

sp.pprint(H.sympyExpr)

b, a = H.getcoeffs()


print(f'\nb coefficients :{b}')
print(f'a coefficients :{a}')

Fs = 48000
f = np.geomspace(1, Fs/2, num=1000)
w = 2*np.pi*f


# Get numerical coefficients for the range of R4 values and C3 values.
# You can also not given a component_values argument to the function.
# In this case, the function will take the default values of the components set in the circuit object.
#3D Case :
'''
component_values = {'C3': np.array([10e-12, 22e-12, 50e-12]),
                    'R4': np.array([4.7e3, 10e3, 22e3, 47e3, 100e3, 500e3])}

b_num, a_num =  H.numerical_analog_coeffs(component_values, combination='parallel')
print(b_num.shape)

h = np.array([[freqs(b_num[i][j], a_num[i][j], worN=w)[1] for j in range(len(a_num[i]))] for i in range(len(b_num))])
plotTransfertFunction(f, h, legend = component_values, semilogx=True, dB=True, phase=True)
'''

#, 'R2': np.array([4.7e3, 10e3, 22e3, 47e3, 100e3, 220e3])
#2D Case
component_values = {'R4': np.array([4.7e3, 10e3, 22e3, 47e3, 100e3, 220e3])}
b_num, a_num =  H.getcoeffs(component_values, z=None, combination='parallel')

h_analog = np.array([freqs(b_num[i], a_num[i], worN=w)[1] for i in range(len(a_num))])
plotTransfertFunction(f, h_analog, legend = component_values, semilogx=True, dB=True, phase=True)


b_num_dig, a_num_dig =  H.getcoeffs(values='num', z='blnr', Fs = Fs, combination='parallel')


#we should not have to convert the coefficients to float64, it should be done automatically...
b_num_dig = np.array(b_num_dig, dtype=np.float64) 
a_num_dig = np.array(a_num_dig, dtype=np.float64)

print(f'\nb coefficients DIGIT :{b_num_dig}')
print(f'a coefficients DIGIT:{a_num_dig}')



#h_discrete = np.array([freqz(b_num_discrete[i], a_num_discrete[i], worN=f)[1] for i in range(len(a_num_discrete))])
_, h_discrete = freqz(b_num_dig, a_num_dig, worN=f, fs=Fs)


fig, ax = plt.subplots(2, 1, figsize=(10, 10))
ax[0].semilogx(f, 20*np.log10(np.abs(h_analog[0])), label='Analog')
ax[0].semilogx(f, 20*np.log10(np.abs(h_discrete)),'--', label='Discrete')
ax[1].semilogx(f, np.angle(h_analog[0]), label='Analog')
ax[1].semilogx(f, np.angle(h_discrete), '--', label='Discrete')

plt.show()


'''
#1D Case
b_num, a_num =  H.numerical_analog_coeffs()

_, h = freqs(b_num, a_num, worN=w)
plotTransfertFunction(f, h, legend='test', semilogx=True, dB=True, phase=True)
'''