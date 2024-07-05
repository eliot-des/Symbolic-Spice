import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

from component_class import AdmittanceComponent, Resistance, Capacitor, Inductance, VoltageSource, ExternalVoltageSource, CurrentSource, IdealOPA
import time 
    
class Circuit:
    """
    Represents a circuit with his associated netlist, which is a collection of electronic components and their connections.
    """

    def __init__(self, inputList):
        """
        Initializes a Circuit object.

        Parameters:
        - inputList (list): A list of strings representing the components in the netlist.

        Returns:
        - None
        """
        self.components         = self.create_component_list(inputList)
        
        self.capacitors         = self.get_components_of_type(self.components, Capacitor)
        self.inductors          = self.get_components_of_type(self.components, Inductance)
        self.voltageSources     = self.get_components_of_type(self.components, VoltageSource)
        self.idealOPAs          = self.get_components_of_type(self.components, IdealOPA)

        self.m = len(self.voltageSources) + len(self.idealOPAs)
        self.n = self.get_nodes_nbr()   # Number of nodes including the ground node

        for component in self.components:
            component.circuit = self


    def create_component_list(self, inputList):
        """
        Creates a list of component objects based on the input list.

        Parameters:
        - inputList (list): A list of strings representing the components in the netlist.

        Returns:
        - components (list): A list of component objects.
        """

        idx = 0
        components = []
        for element in inputList:
            component = self.create_component(element, idx)
            if isinstance(component, VoltageSource) or isinstance(component, IdealOPA):
                idx += 1    
            components.append(component)
        return components


    
    def create_component(self, netlist_line, idx):
        """
        Creates a component object based on the input string. Very simple parser.

        Parameters:
        - netlist_line (str): A string representing a component in the netlist.
        - idx (int): The index of the component in the netlist (used to add a column and row to the initial A matrix).

        Returns:
        - component (Component): A component object.
        """

        symbol, start_node, end_node, value = netlist_line.split()

        if symbol.startswith('R'):
            return Resistance(start_node, end_node, symbol, value)
        
        elif symbol.startswith('C'):
            return Capacitor(start_node, end_node, symbol, value)
        
        elif symbol.startswith('L'):
            return Inductance(start_node, end_node, symbol, value)
    
        elif symbol.startswith('V'):
            if symbol.startswith('Vin'):
                return ExternalVoltageSource(start_node, end_node, symbol, value, idx)
            else:
                return VoltageSource(start_node, end_node, symbol, value, idx)
        elif symbol.startswith('I'):
            return CurrentSource(start_node, end_node, symbol, value)
        
        elif symbol.startswith('O'):
            #value is     : the output node of the OPA
            #start node is: the + terminal of the OPA
            #end node is  : the - terminal of the OPA
            output_node = int(value)
            return IdealOPA(start_node, end_node, output_node, value, idx) 
        else:
            raise ValueError(f"Unknown component symbol: {symbol}")



    def get_components_of_type(self, components, component_type):
        """
        Returns the total number of one type of component from the list of components. 
        The arg component_type must be a subclass of the Component Class
        """
        components_list = [component for component in components if isinstance(component, component_type)]
        return components_list
    
    def get_nodes_nbr(self):
        """
        Returns the total number of unique nodes from the list of components. 
        The ground node is included 
        The set data type is used for avoid to count multiple time the same node (no duplication)
        """
        nodes = set()
        for component in self.components:
            nodes.add(component.start_node)
            nodes.add(component.end_node)
        return len(nodes)



    def stamp_system(self):

        if not hasattr(self, 'A') or not hasattr(self, 'x') or not hasattr(self, 'b'):
            #initialize the A magtrix, x and b vectors with zeros
            self.A = sp.zeros(self.n + self.m, self.n + self.m)
            self.x = sp.zeros(self.n + self.m, 1)
            self.b = sp.zeros(self.n + self.m, 1)

        for component in self.components:
            component.stamp()

        #set the symbolic value at each voltage node in the x 'state' vector:
        for i in range(self.n):
            self.x[i] = sp.symbols(f'v{i}')

        #slice the A matrix, b and x vectors to remove the datum node (ground node)
        self.A = self.A[1:, 1:]
        self.b = self.b[1:, :]
        self.x = self.x[1:, :]

    def print_system(self):
        print('\n\nA matrix:\n')
        sp.pprint(self.A)
        print('\n\nx vector:\n')
        sp.pprint(self.x)
        print('\n\nb vector:\n')
        sp.pprint(self.b)

    def solve_system(self, simplify = False, use_symengine = False):
        """
        Solve the linear system of equations Ax = b

        The solution is stored in the x_solution attribute of the Netlist object.

        Parameters:
        - simplify (bool): If True, simplify the solution. Notes that this can be time consuming.
        - use_symengine (bool): If True, use the symengine library to solve the system to speed up the process.

        Returns:
        - None
        """
        if not hasattr(self, 'A') and not hasattr(self, 'b') and not hasattr(self, 'x'):
            self.stamp_system()


        if use_symengine:
            #check if the symengine library is installed
            try:
                import symengine as se
            except ImportError:
                print('Symengine library is not installed. Please install it to use it.')
                return
            
            A_se = se.Matrix(self.A.tolist())
            b_se = se.Matrix(self.b.tolist())

            x_se = A_se.LUsolve(b_se)
            
            self.x_solution = sp.Matrix(x_se.tolist())
        else:
            self.x_solution = self.A.LUsolve(self.b)

        if simplify:
            self.x_solution = sp.simplify(self.x_solution)


    def get_symbolic_transfert_function(self, output_node, input_node):
        """
        Returns a symbolic transfer function object of the netlist.
        The output_node and input_node are the nodes where the output and input are connected.
        """
        #if there is no x_solution attribute, we need to create one:
        if not hasattr(self, 'x_solution'):
            self.solve_system()

        if output_node > 0:
            output_node_symbol = self.x_solution[output_node - 1]
        else:
            output_node_symbol = 0
        
        if input_node > 0:
            input_node_symbol = self.x_solution[input_node - 1]
        else:
            input_node_symbol = 0
        
        #if there is capacitors or inductors in the circuit, the symbolic transfer function can be factorized:
        if self.capacitors or self.inductors:
            #return the symbolic transfer function into the standard canonical form, where the laplace variable's' is the polynomial variable
            H =  sp.cancel(sp.simplify(output_node_symbol / input_node_symbol), sp.symbols('s'))
        else :
            H = sp.simplify(output_node_symbol / input_node_symbol)

        return CircuitSymbolicTransferFunction(H, self.components)

    def display_components(self, components_list = None):
        """
        Print the circuit's components. By default print all the components.
        The argument components_list can be the voltageSources, currentSources or passiveComponents of the circuit.
        """
        if components_list is None:
            components_list = self.components

        strLine = '------------------------+--------+-------+-------+-----------------+-------------' 
        print(strLine)
        print(f" {'Type':<22} | {'Symbol':^6} | {'Start':^5} | {'End':^5} | {'Admittance':^15} | {'Value':^11} ")
        print(strLine)
        for component in components_list:
            print(f" {str(component):<22} | {str(component.symbol):^6} | {component.start_node:^5} | {component.end_node:^5} | {str(component.admittance):^15} | {component.value:^11} ")
        print(strLine)


def extract_symbolic_analog_filter_coefficients(H, polynomial_variable = sp.symbols('s')):
    """
    Extract the coefficients of the numerator and denominator of a symbolic transfer function.
    """
    num, den = sp.fraction(H)

    num = sp.Poly(num, polynomial_variable)
    den = sp.Poly(den, polynomial_variable)

    return num.all_coeffs(), den.all_coeffs()



class CircuitSymbolicTransferFunction:
    '''
    This class is used to extract the symbolic transfer function of a circuit object, 
    and must contain all the components attributes of the given circuit associated 
    (in order to do the substitutions of the components values in the symbolic transfer function)

    I want to be able to still use the sympy library to manipulate the symbolic transfer function, 
    that's why I override the __getattr__ method to get the attribute of the symbolic transfer function directly
    '''
    def __init__(self, H, components):
        self.sympyExpr = H
        self.components = components
        self.b, self.a = None, None

    def symbolic_analog_filter_coefficients(self):
        num, denum = sp.fraction(self.sympyExpr)
        self.b = sp.Poly(num, sp.symbols('s')).all_coeffs()
        self.a = sp.Poly(denum, sp.symbols('s')).all_coeffs()
        return self.b, self.a

    def numerical_analog_filter_coefficients(self, polynomial_variable=sp.symbols('s')):

        if self.b is None or self.a is None:
            self.b, self.a = self.sympyExpr.extract_symbolic_analog_filter_coefficients()
           
        substitutions = {component.symbol: component.value for component in self.components}

        #convert sympy float object to numpy float object !!!
        num_substituted = np.array([float(coeff.subs(substitutions)) for coeff in self.b])
        den_substituted = np.array([float(coeff.subs(substitutions)) for coeff in self.a])

        #num_substituted = np.array([coeff.subs(substitutions) for coeff in self.b])
        #den_substituted = np.array([coeff.subs(substitutions) for coeff in self.a])
        return num_substituted, den_substituted
    
    def __str__(self):
        return str(self.sympyExpr)
    
    def __repr__(self):
        return repr(self.sympyExpr)
    
    def __getattr__(self, name):
        return getattr(self.sympyExpr, name)
    
    def __call__(self, *args, **kwargs):
        return self.sympyExpr(*args, **kwargs)




def plotTransfertFunction(f, h, title, semilogx=True, dB=True, phase=True):
    ''' 
    Plot the frequency response of the filter

    Parameters:
    - f (array): frequency array (Hz)
    - h (array): frequency response array (complex)
    - title (str): title of the plot
    - semilogx (bool): If True, use a logarithmic scale for the x-axis
    - semilogy (bool): If True, use a logarithmic scale for the y-axis
    - dB (bool): If True, plot the magnitude in dB
    - phase (bool): If True, plot the phase in radians
    '''

    if dB: 
        h_mag =  20 * np.log10(abs(h))
        scale='dB'
    else:
        h_mag = np.abs(h)
        scale='lin.'

    fig, ax = plt.subplots(2 if phase else 1, 1, sharex=True, layout='tight')
    ax = np.atleast_1d(ax)

    ax[0].plot(f, h_mag)

    ax[0].set_title(title)
    ax[0].set_ylabel(f'Magnitude [{scale}]')
    ax[0].grid(which='both')
    ax[-1].set_xticks([1, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000])
    ax[-1].set_xticklabels(['1', '20', '50', '100', '200', '500', '1k', '2k', '5k', '10k', '20k'])
    ax[-1].set_xlim([f[0], f[-1]])
    ax[0].format_coord = lambda x, y: f'x={x:.3f}, y={y:.3f}'

    if semilogx:
        ax[0].set_xscale('log')

    if phase:
        ax[1].plot(f, np.angle(h))
        
        ax[1].set_ylabel('Phase [radians]')
        ax[1].set_xlabel('Frequency [Hz]')
        ax[1].set_yticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
        ax[1].set_yticklabels(['$-\pi$', '$-\pi/2$', '0', '$\pi/2$', '$\pi$'])
        ax[1].set_xlim([f[0], f[-1]])
        ax[1].format_coord = lambda x, y: f'x={x:.3f}, y={y:.2f}'
        ax[1].grid(which='both')

    plt.show()
