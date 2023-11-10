# -*- coding: utf-8 -*-
"""
Created on Sat Oct 21 23:29:53 2023

@author: eliot
"""
import sympy as sp


class Component:
    def __init__(self, start_node, end_node, symbol):
        self.start_node = start_node
        self.end_node = end_node
        self.symbol = sp.symbols(symbol)  # Will be set by child classes
        self.admittance = None
        self.voltage = None

    def __str__(self):
        return type(self).__name__
    
    def get_voltage(self, x_vector):
        
        x_vector = sp.zeros(1,1).col_join(x_vector)
        v_start_node = x_vector[self.start_node]
        v_end_node = x_vector[self.end_node]
        
        voltage = v_end_node - v_start_node
        
        return sp.simplify(voltage)
        

class VoltageSource(Component):
    def __init__(self, start_node, end_node, index):
        super().__init__(start_node, end_node, index)
        
        
class CurrentSource(Component):
    def __init__(self, start_node, end_node, index):
        super().__init__(start_node, end_node, index)


class Resistance(Component):
    def __init__(self, start_node, end_node, index):
        super().__init__(start_node, end_node, index)
        
        R = self.symbol
        self.admittance = 1/R


class Capacitor(Component):
    def __init__(self, start_node, end_node, index):
        super().__init__(start_node, end_node, index)

        omega, C = sp.symbols('omega'), self.symbol
        self.admittance = sp.I * omega * C


class Inductance(Component):
    def __init__(self, start_node, end_node, index):
        super().__init__(start_node, end_node, index)
        
        omega, L = sp.symbols('omega'), self.symbol
        self.admittance = 1 / (sp.I * omega * L)

    
    
    
    