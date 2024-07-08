# -*- coding: utf-8 -*-
"""
Example of using a LT Spice
netlist of the FMV Tone Stack
to retrieve the transfer function
& frequency response.

Created on 07/08/2024
@author: Guillermo A. R.
"""


from scipy.signal import freqs
import numpy as np
import sympy as sp
from symbolicspice import Circuit, plotTransfertFunction
import matplotlib.pyplot as plt

#declare a circuit object
circuit = Circuit(r'FMV Tone Stack.net')

# Get the symbolic transfer function between the output node and the input node
# The input is always set as the first node and the ouput if labeled as the last node
H = circuit.get_symbolic_transfert_function(output_node = 6, input_node = 1)

print('\n\nSymbolic transfer function:\n')
sp.pprint(H.sympyExpr)

# Get symbolic analog coeffs
b, a = H.symbolic_analog_filter_coefficients()

print(f'\nb coefficients :{b}')
print(f'a coefficients :{a}')

# Create freq axis
f = np.geomspace(1,20000, num=1000)
w = 2*np.pi*f


# Get numerical analog coefficients
b_num, a_num =  H.numerical_analog_filter_coefficients()

# Import LT Spice Results
data = open('FMV Tone Stack FR.txt','r')
lines = data.readlines()[1:]
data.close()
data = np.char.replace(lines, '\t', ',')
data = np.char.replace(data, '\n', '')
data = np.char.split(data, ',')

spiceF = np.zeros(len(data))
spiceH = np.zeros(len(data), dtype='complex')

for n, sample in enumerate(data):
    spiceF[n], spiceH[n]= float(sample[0]), float(sample[1]) + 1j * float(sample[2])

_, h = freqs(b_num, a_num, worN=w)

plt.figure(figsize=(6,4))
plt.subplot(2,1,1)
plt.semilogx(spiceF, 20 * np.log10( np.abs(spiceH) ), label='LT Spice')
plt.semilogx(w/2/np.pi, 20 * np.log10( np.abs(h) ), label='Symbolic-Spice', linestyle='--')
plt.xlim((1,20e3))
plt.xlabel('Frequency [Hz]')
plt.ylabel('Magnitude [dB]')
plt.grid()
plt.legend()
plt.subplot(2,1,2)
plt.semilogx(spiceF, np.angle(spiceH), label='LT Spice')
plt.semilogx(w/2/np.pi, np.angle(h), label='Symbolic-Spice', linestyle='--')
plt.xlim((1,20e3))
plt.xlabel('Frequency [Hz]')
plt.ylabel('Phase [rad/s]')
plt.grid()
plt.legend()
plt.suptitle('Frequency Response Comparison')
plt.tight_layout()
plt.show()