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

    def stamp_MNA(self):
        raise NotImplementedError("This method should be implemented by subclasses.")


class AdmittanceComponent(Component):
    def __init__(self, start_node, end_node, symbol, value):
        super().__init__(start_node, end_node, symbol, value)
        self.admittance = None # This is set in subclasses

    def stamp_MNA(self):
        #method to stamp the admittance of the components in the A matrix of the circuit class
        self.circuit.A[self.start_node, self.start_node] += self.admittance
        self.circuit.A[self.end_node,     self.end_node] += self.admittance    
        # Off-diagonal elements
        self.circuit.A[self.start_node,   self.end_node] -= self.admittance
        self.circuit.A[self.end_node,   self.start_node] -= self.admittance
    
    def stamp_MNA_ss(self, A, x, b, DLC, SLC, SIV0):
        #method to stamp the admittance of the components in a MNA system 
        #that will be used to derive the state-space representation
        A[self.start_node, self.start_node] += self.admittance
        A[self.end_node,     self.end_node] += self.admittance
        A[self.start_node,   self.end_node] -= self.admittance
        A[self.end_node,   self.start_node] -= self.admittance
        

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

        #only needed for the state-space representation
        self.voltage = sp.symbols(f'v_{self.symbol}')
        self.current = sp.symbols(f'i_{self.symbol}')

    def index_ss(self, index):
        self.index_ss = index
        
    def stamp_MNA_ss(self, A, x, b, DLC, SLC, SIV0):
        #method to stamp the capacitor in a MNA system used to derive the state-space representation
        #Here the capacitor is represented as a pure voltage source

        n = self.circuit.n
        nC = self.circuit.nC
        nL = self.circuit.nL
        
        A[n + self.index_ss, self.start_node] =  1
        A[n + self.index_ss,   self.end_node] = -1
        A[self.start_node, n + self.index_ss] =  1
        A[self.end_node,   n + self.index_ss] = -1

        # Stamp the b vector for this voltage source
        b[n + self.index_ss] = self.voltage

        # Stamp the unknown current accross the voltage source in the x 'state' vector
        x[n + self.index_ss] = self.current

        SLC[n + self.index_ss, nC+nL] = 1
        DLC[nC+nL, n + self.index_ss] = 1/self.symbol

        self.circuit.nC+=1
    

class Inductance(AdmittanceComponent):
    def __init__(self, start_node, end_node, symbol, value):
        super().__init__(start_node, end_node, symbol, value)
        s, L = sp.symbols('s'), self.symbol
        self.admittance = 1 / (s * L)


        #only needed for the state-space representation
        self.voltage = sp.symbols(f'v_{self.symbol}')
        self.current = sp.symbols(f'i_{self.symbol}')

    def stamp_MNA_ss(self, A, x, b, DLC, SLC, SIV0):
        #method to stamp the inductor in a MNA system used to derive the state-space representation
        #Here the inductor is represented as a pure current source
        n = self.circuit.n
        nC = self.circuit.nC
        nL = self.circuit.nL
        
        b[self.start_node] -= self.current
        b[self.end_node  ] += self.current

        SLC[self.start_node, nC+nL] = -1
        SLC[self.end_node, nC+nL] = 1

        DLC[nC+nL, self.start_node] = 1/self.symbol
        DLC[nC+nL, self.end_node] = -1/self.symbol

        self.circuit.nL+=1
        




class VoltageSource(Component):
    def __init__(self, start_node, end_node, symbol, value, index):
        super().__init__(start_node, end_node, symbol, value)
        self.index = index
        
        self.voltage = self.symbol
        self.current = sp.symbols(f'i_{self.symbol}')

    def stamp_MNA(self):
        n = self.circuit.n
        
        self.circuit.A[n + self.index, self.start_node] =  1
        self.circuit.A[n + self.index,   self.end_node] = -1
        self.circuit.A[self.start_node, n + self.index] =  1
        self.circuit.A[self.end_node,   n + self.index] = -1

        # Stamp the b vector for this voltage source
        self.circuit.b[n + self.index] = self.voltage

        # Stamp the unknown current accross the voltage source in the x 'state' vector
        self.circuit.x[n + self.index] = self.current
    
    def index_ss(self, index):
        self.index_ss = index

    def stamp_MNA_ss(self, A, x, b, DLC, SLC, SIV0):
        n = self.circuit.n
        
        A[n + self.index_ss, self.start_node] =  1
        A[n + self.index_ss,   self.end_node] = -1
        A[self.start_node, n + self.index_ss] =  1
        A[self.end_node,   n + self.index_ss] = -1

        # Stamp the b vector for this voltage source
        b[n + self.index_ss] = self.symbol

        # Stamp the unknown current accross the voltage source in the x 'state' vector
        x[n + self.index_ss] = self.current


class ExternalVoltageSource(VoltageSource):
    def __init__(self, start_node, end_node, symbol, value, index):
        super().__init__(start_node, end_node, symbol, value, index)

    def stamp_MNA_ss(self, A, x, b, DLC, SLC, SIV0):
        #call the stamp_MNA_ss method of the parent class
        super().stamp_MNA_ss(A, x, b, DLC, SLC, SIV0)

        nI = self.circuit.nI
        nV = self.circuit.nV

        SIV0[n, nI+nV] = 1
        nV+=1

class VoltageProbe(Component):

    def __init__(self, start_node, end_node, symbol, value, index):
        super().__init__(start_node, end_node, symbol, value)
        self.index = index



class CurrentSource(Component):
    def __init__(self, start_node, end_node, symbol, value):
        super().__init__(start_node, end_node, symbol, value)
        self.current = self.symbol

    def stamp_MNA(self):
        self.circuit.b[self.start_node] -= self.current
        self.circuit.b[self.end_node  ] += self.current
    
    def stamp_MNA_ss(self, A, x, b, DLC, SLC, SIV0):
        b[self.start_node] -= self.current
        b[self.end_node  ] += self.current

class ExternalCurrentSource(CurrentSource):
    def __init__(self, start_node, end_node, symbol, value):
        super().__init__(start_node, end_node, symbol, value)

    def stamp_MNA_ss(self, A, x, b, DLC, SLC, SIV0):
        super().stamp_MNA_ss(A, x, b, DLC, SLC, SIV0)

        nI = self.circuit.nI
        nV = self.circuit.nV

        SIV0[self.start_node, nI+nV] = -1
        SIV0[self.end_node, nI+nV] = 1

        nI+=1

class IdealOPA(Component):
    def __init__(self, start_node, end_node, output_node, symbol, index):
        super().__init__(start_node, end_node, symbol)

        #start node is: the + terminal of the OPA
        #end node is  : the - terminal of the OPA
        self.output_node = int(output_node)
        self.index = index

    def stamp_MNA(self):
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

    def stamp_MNA(self):
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

    def stamp_MNA(self):
        n = self.circuit.n

        self.circuit.A[self.start_node, self.start_node_2] = self.symbol
        self.circuit.A[self.start_node,   self.end_node_2] = -self.symbol
        self.circuit.A[self.end_node,   self.start_node_2] = -self.symbol
        self.circuit.A[self.end_node,     self.end_node_2] = self.symbol

        self.circuit.A[self.start_node_2, self.start_node] = -self.symbol
        self.circuit.A[self.start_node_2,   self.end_node] = self.symbol
        self.circuit.A[self.end_node_2,   self.start_node] = self.symbol
        self.circuit.A[self.end_node_2,     self.end_node] = -self.symbol