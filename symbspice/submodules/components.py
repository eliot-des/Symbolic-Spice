# -*- coding: utf-8 -*-
"""
Created on Thur Jul 04 17:29:58 2024

@author: eliot

The value of the reactive components are set by linearizing the 
components thanks to the Backward Euler Method
"""

import sympy as sp

class Component:
    def __init__(self, start_node, end_node, symbol, value=0):
        self.start_node = int(start_node)
        self.end_node = int(end_node)
        self.symbol = sp.symbols(symbol)
        self.value = float(value)
        self.admittance = None  # This will be set in subclass

    def __str__(self):
        return type(self).__name__
    
    def setup(self, circuit):
        self.circuit = circuit

    def stamp_MNA():
        raise NotImplementedError("This method should be implemented by subclasses.")


class AdmittanceComponent(Component):
    def __init__(self, start_node, end_node, symbol, value):
        super().__init__(start_node, end_node, symbol, value)
        self.admittance = None # This is set in subclasses

    def stamp_MNA():
        #method to stamp the admittance of the components in the A matrix of the circuit class
        self.circuit.A[self.start_node, self.start_node] += self.admittance
        self.circuit.A[self.end_node,     self.end_node] += self.admittance    
        # Off-diagonal elements
        self.circuit.A[self.start_node,   self.end_node] -= self.admittance
        self.circuit.A[self.end_node,   self.start_node] -= self.admittance

class Resistance(AdmittanceComponent):     
    def __init__(self, start_node, end_node, symbol, value):
        super().__init__(start_node, end_node, symbol, value)
        R = self.symbol
        self.admittance = 1 / R


class Capacitor(AdmittanceComponent):
    def __init__(self, start_node, end_node, symbol, value):
        super().__init__(start_node, end_node, symbol, value)
        s, C = sp.symbols('s'), self.symbol
        self.admittance = s * C
    

class Inductance(AdmittanceComponent):
    def __init__(self, start_node, end_node, symbol, value):
        super().__init__(start_node, end_node, symbol, value)
        s, L = sp.symbols('s'), self.symbol
        self.admittance = 1 / (s * L)



class VoltageSource(Component):
    def __init__(self, start_node, end_node, symbol, value, index):
        super().__init__(start_node, end_node, symbol, value)
        self.index = index

    def stamp_MNA():
        n = self.circuit.n
        
        self.circuit.A[n + self.index, self.start_node] =  1
        self.circuit.A[n + self.index,   self.end_node] = -1
        self.circuit.A[self.start_node, n + self.index] =  1
        self.circuit.A[self.end_node,   n + self.index] = -1

        # Stamp the b vector for this voltage source
        self.circuit.b[n + self.index] = self.symbol

        # Stamp the unknown current accross the voltage source in the x 'state' vector
        self.circuit.x[n + self.index] = sp.symbols(f'i_{self.symbol}')


class ExternalVoltageSource(VoltageSource):
    def __init__(self, start_node, end_node, symbol, value, index):
        
        super().__init__(start_node, end_node, symbol, value, index)

class VoltageProbe(Component):

    def __init__(self, start_node, end_node, symbol, value, index):
        super().__init__(start_node, end_node, symbol, value)
        self.index = index



class CurrentSource(Component):
    def __init__(self, start_node, end_node, symbol, value):
        super().__init__(start_node, end_node, symbol, value)

    def stamp_MNA():
        self.circuit.b[self.start_node] -= self.symbol
        self.circuit.b[self.end_node  ] += self.symbol

class IdealOPA(Component):
    def __init__(self, start_node, end_node, output_node, symbol, index):
        super().__init__(start_node, end_node, symbol)

        #start node is: the + terminal of the OPA
        #end node is  : the - terminal of the OPA
        self.output_node = int(output_node)
        self.index = index

    def stamp_MNA():
        n = self.circuit.n 
        
        self.circuit.A[self.output_node, n + self.index] =  1
        self.circuit.A[n + self.index,  self.start_node] =  1
        self.circuit.A[n + self.index,    self.end_node] = -1

class Transformer(Component):
    def __init__(self, start_node, end_node, start_node_2, end_node_2, symbol, value):
        super().__init__(start_node, end_node, symbol, value)
        self.start_node_2 = start_node_2
        self.end_node_2 = end_node_2
        #here, the value of the transformer is the wires ratio

    def stamp_MNA():
        n = self.circuit.n

        self.circuit.A[self.start_node,   n + self.index] =  1
        self.circuit.A[self.end_node,     n + self.index] = -1
        self.circuit.A[self.start_node_2, n + self.index] = self.symbol
        self.circuit.A[self.end_node_2,   n + self.index] = -self.symbol

        self.circuit.A[n + self.index,   self.start_node] =  1
        self.circuit.A[n + self.index,     self.end_node] = -1
        self.circuit.A[n + self.index, self.start_node_2] = -self.symbol
        self.circuit.A[n + self.index,   self.end_node_2] = self.symbol

        self.circuit.b[n + self.index] = sp.symbols(f'i_{self.symbol}')

class Gyrator(Component):
    def __init__(self, start_node, end_node, start_node_2, end_node_2, symbol, value):
        super().__init__(start_node, end_node, symbol, value)
        self.start_node_2 = start_node_2
        self.end_node_2 = end_node_2
        #here, the value of the gyrator is the transconductance

    def stamp_MNA():
        n = self.circuit.n

        self.circuit.A[self.start_node, self.start_node_2] = self.symbol
        self.circuit.A[self.start_node,   self.end_node_2] = -self.symbol
        self.circuit.A[self.end_node,   self.start_node_2] = -self.symbol
        self.circuit.A[self.end_node,     self.end_node_2] = self.symbol

        self.circuit.A[self.start_node_2, self.start_node] = -self.symbol
        self.circuit.A[self.start_node_2,   self.end_node] = self.symbol
        self.circuit.A[self.end_node_2,   self.start_node] = self.symbol
        self.circuit.A[self.end_node_2,     self.end_node] = -self.symbol