# -*- coding: utf-8 -*-
"""
Created on Sat Oct 21 23:29:53 2023

@author: eliot
"""
import sympy as sp


class Component:
    def __init__(self, start_node, end_node, symbol, value):
        self.start_node = start_node
        self.end_node = end_node
        self.symbol = sp.symbols(symbol)
        self.value = value
        self.admittance = None  # This will be set in each subclass

    def __str__(self):
        return type(self).__name__

class PassiveComponent(Component):
    def __init__(self, start_node, end_node, symbol, value):
        super().__init__(start_node, end_node, symbol, value)
        

class Resistance(PassiveComponent):
    def __init__(self, start_node, end_node, index, value):
        super().__init__(start_node, end_node, index, value)
        R = self.symbol
        self.admittance = 1 / R

class Capacitor(PassiveComponent):
    def __init__(self, start_node, end_node, index, value):
        super().__init__(start_node, end_node, index, value)
        s, C = sp.symbols('s'), self.symbol
        self.admittance = s * C

class Inductance(PassiveComponent):
    def __init__(self, start_node, end_node, index, value):
        super().__init__(start_node, end_node, index, value)
        s, L = sp.symbols('s'), self.symbol
        self.admittance = 1 / (s * L)

class VoltageSource(Component):
    def __init__(self, start_node, end_node, index, value):
        super().__init__(start_node, end_node, index, value)

class CurrentSource(Component):
    def __init__(self, start_node, end_node, index, value):
        super().__init__(start_node, end_node, index, value)

    
    
    
    
