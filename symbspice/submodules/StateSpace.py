import sympy as sp
import numpy as np



class StateSpace:
    def __init__(self, circuit, A, x, B, u, C, y, D):
        self.circuit = circuit
        self.A = A
        self.x = x
        self.B = B
        self.u = u
        self.C = C
        self.y = y
        self.D = D

    def state_equations(self):
        """
        Show the state equation of the system in the console

        dx/dt = Ax + Bu
        """
        
        if len(self.x) == 0:
            return print('The circuit has no state variables...')
        
        x_dot = sp.Matrix([sp.Derivative(xi, 't') for xi in self.x])

        if len(self.u) == 0:
            return sp.Eq(x_dot, sp.MatMul(self.A, self.x))
        else:
            return sp.Eq(x_dot, sp.MatAdd(sp.MatMul(self.A, self.x) ,sp.MatMul(self.B, self.u)))

    def output_equations(self):
        """
        Show the output equation of the system in the console

        y = Cx + Du
        """
        
        if len(self.y) == 0:
            return print('The circuit has no output variables...')
        
        if len(self.u) == 0:
            return sp.Eq(self.y, sp.MatMul(self.C, self.x))
        else:
            return sp.Eq(self.y, sp.MatAdd(sp.MatMul(self.C, self.x), sp.MatMul(self.D, self.u)))
        