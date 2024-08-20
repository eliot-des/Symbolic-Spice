import sympy as sp
import numpy as np



class StateSpace:
    def __init__(self, circuit, A, x, B, u, C, D, y):
        self.circuit = circuit
        self.A = A
        self.x = x
        self.B = B
        self.u = u
        self.C = C
        self.D = D
        self.y = y