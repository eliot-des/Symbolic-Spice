import sympy as sp
import numpy as np
import itertools
import matplotlib.pyplot as plt
from numbers import Number
import warnings

from components import AdmittanceComponent, Resistance, Capacitor, Inductance, VoltageSource, ExternalVoltageSource, CurrentSource, IdealOPA, Transformer, Gyrator
import time 
    
class Circuit:
    """
    Represents a circuit with his associated netlist, which is a collection of electronic components and their connections.
    """

    def __init__(self, inputList):
        """
        Initializes a Circuit object.

        Parameters:
        - inputList (list or string): A list of strings representing the components 
        in the netlist or a path to a LT Spice netlist file.

        Returns:
        - None
        """

        # If input is netlist file, parse netlist
        if isinstance(inputList, str):
            inputList = self.__loadnet(inputList)
        
        self.components         = self.__create_component_list(inputList)
        
        self.capacitors         = self.__get_components_of_type(self.components, Capacitor)
        self.inductors          = self.__get_components_of_type(self.components, Inductance)
        self.voltageSources     = self.__get_components_of_type(self.components, VoltageSource)
        self.idealOPAs          = self.__get_components_of_type(self.components, IdealOPA)

        self.m = len(self.voltageSources) + len(self.idealOPAs)
        self.n = self.get_nodes_nbr()   # Number of nodes including the ground node

        for component in self.components:
            component.circuit = self


    def __create_component_list(self, inputList):
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
            component = self.__create_component(element, idx)
            if isinstance(component, VoltageSource) or isinstance(component, IdealOPA):
                idx += 1    
            components.append(component)
        return components


    
    def __create_component(self, netlist_line, idx):
        """
        Creates a component object based on the input string. Very simple parser.

        Parameters
        ----------
        - netlist_line (str): A string representing a component in the netlist.
        - idx (int): The index of the component in the netlist (used to add a column and row to the initial A matrix).

        Returns
        -------
        - component (Component): A component object.
        """

        component_map = {
            'R': Resistance,
            'C': Capacitor,
            'L': Inductance,
            'V': VoltageSource,
            'I': CurrentSource,
            'O': IdealOPA,
            'T': Transformer,
            'G': Gyrator
        }

        # Split the netlist line into its components
        parts = netlist_line.split()
        symbol = parts[0]
        args = parts[1:]

        # Determine the type of component and instantiate it. Not a big deel to do a for loop 
        # and check if the symbol starts with the key of the component_map because the number 
        # of components is very limited.
        for key, component_class in component_map.items():
            if symbol.startswith(key):
                if key == 'V':
                    if symbol.startswith('Vin'):
                        return ExternalVoltageSource(*args[:-1], symbol, args[-1], idx)
                    else:
                        return component_class(*args[:-1], symbol, args[-1], idx) #normal voltage source
                if key == 'O':
                    return component_class(*args, symbol, idx)
                return component_class(*args[:-1], symbol, args[-1])
                
        raise ValueError(f"Unknown component symbol: {symbol}")



    def __get_components_of_type(self, components, component_type):
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
        #initialize the A matrix, x and b vectors with zeros
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

        Parameters
        ----------
        - simplify (bool): If True, simplify the solution. Notes that this can be time consuming.
        - use_symengine (bool): If True, use the symengine library to solve the system to speed up the process.

        Returns
        -------
        - None
        """
        if not hasattr(self, 'A') or not hasattr(self, 'b') or not hasattr(self, 'x'):
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


    def tf(self, output_node, input_node=1, norm = True):
        """
        Returns the b & a coefficients of the symbolic transfer function object of the netlist.

        Parameters
        ----------
        output_node (int): Node where the output is set.
        input_node (int): Node where Vin is set.
        norm (bool): Normalize coefficients (b/a0 & a/a0 if True).
        """
        # Ensure the system is solved:
        if not hasattr(self, 'x_solution'):
            self.solve_system()

        # Function to retrieve node symbol:
        def get_Voltage_Expr(node):
            if node == 0:
                return 0
            elif node <= self.n:
                return self.x_solution[node - 1]
            else:
                raise ValueError(f"Node {node} is not in the circuit.")

        # Retrieve symbols for output and input nodes:
        output_Expr = get_Voltage_Expr(output_node)
        input_Expr  = get_Voltage_Expr(input_node)

        # Calculate the symbolic transfer function:
        H = output_Expr / input_Expr

        if self.capacitors or self.inductors:
            H = sp.cancel(sp.simplify(H), sp.symbols('s'))
        else:
            H = sp.simplify(H)


        transfer_function = TransferFunction(H, self.components)

        if norm: transfer_function.normalized()

        return transfer_function


    def __loadnet(self, fname):
        """
        Converts a LT Spice netlist into a list
        that follows the input format of the Circuit class.

        Parameters:
        - fname (String): Path to circuit's netlist file.

        Returns:
        - outlist (list): list of formatted netlist components.
        """ 
        # Netlist import memory
        temp = []
        with open(fname) as text:
            for line in text:
                # purge lines after netlist
                if line.startswith('.'):
                    break

                # Modify node names to match parsing
                line = line.replace('N00', 'N0')
                line = line.replace('N0', '')

                # Remove special/unnecesary characters
                line = line.replace(' AC', '')
                
                # separate by lines
                temp.append(line.strip('\n'))

        # remove header
        temp.pop(0)

        # Create proper netlist matrix
        netlist = np.zeros((len(temp), 4), dtype='<U22')
        # set values to zero
        netlist[:,3] = str(0)
        # Separate name, start node, end node, value
        for count, component in enumerate(temp):
            # if value is declared
            if np.shape(component.split(' '))[0] == 4:
                netlist[count,:] = component.split(' ')
            # else, value is set to 0
            else:
                netlist[count,:3] = component.split(' ')

        # Replace SPICE units with Python units
        netlist[:,3] = np.char.replace(netlist[:,3], 'p', 'e-12')
        netlist[:,3] = np.char.replace(netlist[:,3], 'n', 'e-9')
        netlist[:,3] = np.char.replace(netlist[:,3], 'u', 'e-6')
        netlist[:,3] = np.char.replace(netlist[:,3], 'μ', 'e-6')
        netlist[:,3] = np.char.replace(netlist[:,3], 'm', 'e-3')
        netlist[:,3] = np.char.replace(netlist[:,3], 'k', 'e3')
        netlist[:,3] = np.char.replace(netlist[:,3], 'MEG', 'e6')

        # Remove unnecessary units
        netlist[:,3] = np.char.replace(netlist[:,3], 'V', '')

        # Find nodes with numbers
        n1 = np.char.isdigit(netlist[:,1:3])
        # Add 1 to max node (in case output is labeled)
        if np.any(netlist[:, 1:3] == 'out'):
            self.outNode = np.max(netlist[:,1:3][n1].astype(int) + 1)
        else: # if not choose last node as output
            self.outNode = np.max(netlist[:,1:3][n1].astype(int))

        # Make sure there is a voltage source from in to ground called Vin
        if np.any(netlist[:, 0] == 'Vin'):
            # Make sure it's value is 1
            node = np.argwhere(netlist[:,0] == 'Vin')[0]
            netlist[node, 3] = '1'
            # Set input to the first node
            netlist[:,1:3] = np.char.replace(netlist[:,1:3], '1', 'temp')
            netlist[:,1:3] = np.char.replace(netlist[:,1:3], netlist[node,1], '1')
            netlist[:,1:3] = np.char.replace(netlist[:,1:3], 'temp', netlist[node,1])
            # Set Vin at the first row
            temp = netlist[0,:].copy()
            netlist[0,:] = netlist[node].flatten()
            netlist[int(node), :] = temp
            # Sort elements by node 1
            netlist = netlist[netlist[:, 1].argsort()]
        else:
            raise Exception("Netlist is missing Vin, add a voltage source named Vin where the input is wanted")

        # If there is a net label for the output set it as the last node
        if np.any(netlist[:, 1:3] == 'out'):
            # Replace out to max node
            netlist[:,1:3] = np.char.replace(netlist[:,1:3], 'out', str(self.outNode))

        # return netlist as list of rows (each row is a component: [name, start node, end node, value])
        outlist = []
        for n in range(netlist.shape[0]):
            # Clean non-numeric values such as R, C, etc
            if np.char.isdigit(testlist[n,3][0]) == False:
                netlist[n,3] = '0'
                
            outlist.append(' '.join(str(item) for item in netlist[n,:]))

        return outlist



class TransferFunction:
    '''
    This class is used when using the get_symbolic_transfert_function() method of the Circuit class,
    in order to extract the symbolic transfer function of a circuit object, between two nodes. 
    The object created must contain all the components attributes of the given circuit associated 
    (in order to do the substitutions of the components values in the symbolic transfer function).

    Normally, you should not need to use this class directly, but only through the 
    get_symbolic_transfert_function() method of the Circuit class, that will return an object of 
    this class.
    '''

    def __init__(self, H, components):
        self.sympyExpr = H
        self.components = components
        self.b, self.a = None, None  #Symbolic analog filter coefficients

    def sym_coeffs(self):
        if self.b is None or self.a is None:
            num, denum = sp.fraction(self.sympyExpr)
            self.b = sp.Poly(num, sp.symbols('s')).all_coeffs()
            self.a = sp.Poly(denum, sp.symbols('s')).all_coeffs()
        return self.b, self.a

    def normalized(self):
        '''
        Normalize the transfer function to have self.a[0] = 1.
        '''
        num, den = sp.fraction(self.sympyExpr)

        if self.b is None or self.a is None:
            self.b, self.a = self.sym_coeffs()

        a0 = self.a[0]
        self.a = [coeff/a0 for coeff in self.a]
        self.b = [coeff/a0 for coeff in self.b]

        self.sympyExpr = sp.Poly(self.b, sp.symbols('s')) / sp.Poly(self.a, sp.symbols('s'))

    def num_coeffs(self, component_values=None, combination='nested'):
        """
        Return the numerical coefficients `b_num` and `a_num` of the analog filter transfer function.
        The coefficients are calculated by substituting the component values in the symbolic transfer function.

        Parameters
        ----------
        component_values : {None, dict}, optional
            A dictionary of component values. The keys are the component symbols, and the values are the component values.
            * If component_values is None, the default values of the components set in the Circuit object will be used.
            
            * If component_values is a dictionary with one key and a 1D array key values, 
            the function will return the a and b numerical coefficients for each value of the array, in a 2D array.
            
            * If component_values is a dictionary with multiple keys and with their associated 1D array key values,
            the function will return the a and b numerical coefficients for each combination of values in a (N+1)-D array,
            where N is the number of keys in the dictionary.

        combination : {'nested', 'parallel'}, optional
            * If the `component_values` dictionary has multiple keys, the 'combinations' argument specifies how to combine the values.
            * If the `component_values` dictionary has only one key, the 'combinations' argument is ignored 
            (for the moment not really, but it should works with 'nested' or 'parallel' for one key too, even if it's not optimal)

            * If `'nested'`, return all the possible combinations of the component values. 
            For example, if the dictionary is `{'R1': [1, 2], 'R2': [3, 4, 5]}`, 
            the `'nested'` option will return the combinations `[[(1, 3), (1, 4), (1, 5)], [(2, 3), (2, 4), (2,5)]]`.
            * If 'parallel', return only the combinations of the component values in the order of the dictionary keys.
            Therefore, the same number of values must be provided for each key in the dictionary.
            For example, if the dictionary is `{'R1': [1, 2], 'R2': [3, 4]}`, 
            the 'parallel' option will return the combinations `[(1, 3), (2, 4)]`.

        Returns
        -------
        b_num : list
            The numerator coefficients of the transfer function.
        a_num : list
            The denominator coefficients of the transfer function.
        """

        if self.b is None or self.a is None:
            self.b, self.a = self.sym_coeffs()
        
        if component_values is None:
            component_values = {component.symbol: component.value for component in self.components}

            b_num = np.array([float(coeff.subs(component_values)) for coeff in self.b])
            a_num = np.array([float(coeff.subs(component_values)) for coeff in self.a])

            return b_num, a_num

        keys, values = zip(*component_values.items())

        if combination == 'nested':
            combinations = list(itertools.product(*values))
            shape = [len(values_set) for values_set in values] + [len(self.b)]
        elif combination == 'parallel':
            if not all(len(value) == len(values[0]) for value in values):
                raise ValueError("All values in the component_values dictionary must have the same length if 'parallel' is selected.")
            combinations = list(zip(*values))
            shape = (len(combinations), len(self.b))
        else:
            raise ValueError("The 'combinations' argument must be either 'nested' or 'parallel'.")

        b_num, a_num = [], []
        for comb in combinations:
            substitutions = { component.symbol: comb[keys.index(str(component.symbol))]
                            if str(component.symbol) in keys else component.value for component in self.components}

            b_num_temp = [float(coeff.subs(substitutions)) for coeff in self.b]
            a_num_temp = [float(coeff.subs(substitutions)) for coeff in self.a]

            b_num.append(b_num_temp)
            a_num.append(a_num_temp)

        return np.array(b_num).reshape(shape), np.array(a_num).reshape(shape)


    def getcoeffs(self, values = 'symb', z = None, Fs = None, combination = 'nested'):
        """
        Return b & a coefficients of the transfer function.

        Parameters
        ----------
        - values : 'symb', 'num', or dict like {'R1': [1e3, 2e6], 'R2': [3e2, 4e2]}
            * If  'sym', return the symbolic coefficients.
            * If  'num', return the numerical coefficients for the default values of the components set in the Circuit object.
            * If a dictionary is provided, return the numerical coefficients for the given values.
        - z : {'frwrd', 'bckwrd', 'blnr'}, or None, optional
            The desired discretization scheme.
            * If None, the function will return the analog filter coefficients.
        - Fs : float, optional
            The sampling frequency.
        - combination : {'nested', 'parallel'}, optional
            * If the `component_values` dictionary has multiple keys, the 'combinations' argument specifies how to combine the values.
            * If the `component_values` dictionary has only one key, the 'combinations' argument is ignored
        Returns
        -------
        b, a : np.array
            The b & a coefficients of the transfer function.
        """
        
        # First, we derive the analog filter coefficients:

        if values == 'symb':
            b, a = self.sym_coeffs()
        elif values == 'num' or isinstance(values, dict):
            component_values = None if values == 'num' else values
            b, a = self.num_coeffs(component_values=component_values, combination=combination)
        else:
            raise ValueError("The 'values' argument must be either 'symb', 'num', or a dictionary of component values.")

        
        # Then we use the eval() function to derive the digital filter coefficients (if the user want) 
        # corresponding to the analog filter coefficients, either they are symbolic or numerical.

        if z is not None:
            if values == 'symb':
                if isinstance(Fs, sp.Symbol):
                    b, a = self.coeffz(b, a, z, Fs)
                elif Fs is None:
                    b, a = self.coeffz(b, a, z)
                else:
                    raise ValueError("The sampling frequency must be a symbolic variable if values = 'symb'.")
            else:
                if Fs == None:
                    raise ValueError("Samplerate was not given.")
                b, a = self.coeffz(b, a, z, Fs)
        return b, a

    def gen_coeffz(self, N, scheme='blnr'):
        """
        Returns the symbolic discretized coefficients for a given order N & discretization scheme.

        Parameters:
        - N (int): Order of transfer function.
        - scheme (str): Desired discretization scheme ('frwrd', 'bckwrd', 'blnr').

        Returns:
        - Bd (np.array): Discretized numerator coefficients.
        - Ad (np.array): Discretized denominator coefficients.
        """
        s, z, Ts = sp.symbols('s z T_s')
        b = sp.symbols(f'a_0:{N + 1}')[::-1]
        a = sp.symbols(f'b_0:{N + 1}')[::-1]

        #create a dictionary of the different discretization schemes
        #therefore, additonal schemes can be added easily
        schemes = {
            'frwrd': (z - 1) / Ts,
            'bckwrd': (z - 1) / (Ts * z),
            'blnr': 2 / Ts * (z - 1) / (z + 1)
        }

        if scheme not in schemes:
            raise ValueError("Invalid scheme given")

        ds = schemes[scheme]
        B = sp.Poly(a, s).as_expr().subs(s, ds) * (z + 1)**N
        A = sp.Poly(b, s).as_expr().subs(s, ds) * (z + 1)**N

        H = sp.simplify(B / A)
        B, A = sp.fraction(H)

        B = sp.collect(sp.expand(B / z**N), 1 / z)
        A = sp.collect(sp.expand(A / z**N), 1 / z)

        Bd = [B.coeff(z, -n) for n in range(N + 1)]
        Ad = [A.coeff(z, -n) for n in range(N + 1)]

        return Bd, Ad



    def sub_coeffz(self, b, a, scheme='blnr', srate=sp.Symbol('F_s')):
        """
        Returns the symbolic or real discretized coefficients for a given array of analog coefficients.

        Parameters:
        - b (list/array): Continuous b coefficients.
        - a (list/array): Continuous a coefficients.
        - scheme (str): Discretization scheme ('frwrd', 'bckwrd', 'blnr').
        - srate (float or sympy.Symbol): Samplerate.

        Returns:
        - Bd (np.array): Discretized numerator coefficients.
        - Ad (np.array): Discretized denominator coefficients.
        """
        N = len(b)
        B, A = self.gen_coeffz(N - 1, scheme)
        subs_dict = {'T_s': 1. / srate}

        b, a = b[::-1], a[::-1]
        for n in range(N):
            subs_dict[f'b_{n}'] = b[n]
            subs_dict[f'a_{n}'] = a[n]

        if all(isinstance(coeff, (int, float)) for coeff in b + a):
            out_type=np.float64
        else:
            out_type = object

        Bd = np.array([B[n].subs(subs_dict) for n in range(N)], dtype=out_type)
        Ad = np.array([A[n].subs(subs_dict) for n in range(N)], dtype=out_type)

        return Bd, Ad



    def coeffz(self, b, a, scheme='blnr', srate=sp.Symbol('F_s')):
        """
        Returns the symbolic or real discretized coefficients for a given array of analog coefficients.

        Parameters:
        - b (list/array): Continuous b coefficients.
        - a (list/array): Continuous a coefficients.
        - scheme (str): Discretization scheme ('frwrd', 'bckwrd', 'blnr').
        - srate (float or sympy.Symbol): Samplerate.

        Returns:
        - Bd (np.array): Discretized numerator coefficients.
        - Ad (np.array): Discretized denominator coefficients.
        """
        b = np.asarray(b)
        a = np.asarray(a)


        if b.shape != a.shape:
            raise ValueError("Shapes of b and a must be the same.")

        Bd = np.empty_like(b)
        Ad = np.empty_like(a)

        #we must iter on b and a by keeping b and a as array, because this a list of coefficient associated to a given filter.
        it = np.nditer(b[..., 0], flags=['multi_index'])

        for _ in it:
            idx = it.multi_index
            Bd[idx], Ad[idx] = self.sub_coeffz(b[idx], a[idx], scheme, srate)
            
        return Bd, Ad

    def __str__(self):
        return str(self.sympyExpr)
    
    def __repr__(self):
        return repr(self.sympyExpr)
    
    def __getattr__(self, name):
        '''
        I want to be able to still use the sympy library to manipulate the symbolic transfer function, 
        that's why I override the __getattr__ method to get the attribute of the symbolic transfer function directly
        '''
        return getattr(self.sympyExpr, name)
    
    def __call__(self, *args, **kwargs):
        return self.sympyExpr(*args, **kwargs)



def plotTransfertFunction(f, h, legend=None, title=None, semilogx=True, dB=True, phase=True):
    ''' 
    Plot the frequency response of the filter.

    Parameters:
    - f (array): frequency array (Hz) -> 1D array
    - h (array): frequency response array (complex) -> Max 3D array!
    - legend (dict or str): Legend of the different lines.
                            Use the same argument as component_values in the 
                            num_coeffs() method if 
                            you want a fast legend creation!
    - title (str): title of the plot
    - semilogx (bool): If True, use a logarithmic scale for the x-axis
    - semilogy (bool): If True, use a logarithmic scale for the y-axis
    - dB (bool): If True, plot the magnitude in dB
    - phase (bool): If True, plot the phase in radians
    '''
    
    if h.ndim > 3:
        raise ValueError("The input array h must have at most 3 dimensions.")

    #default colors and linestyles 
    colors = plt.cm.tab10(np.linspace(0, 1, 10)) #RdYlBu, plasma, viridis, inferno, magma, cividis
    linestyles = ['-', '--', '-.', ':', '-*', '-o', '-+', '-x', '-|']
    colors_nbr = len(colors)
    linestyles_nbr = len(linestyles)

    h = np.atleast_2d(h)

    if dB: 
        h_mag =  20 * np.log10(np.abs(h))
        scale = 'dB'
    else:
        h_mag = np.abs(h)
        scale = 'lin.'

    fig, ax = plt.subplots(2 if phase else 1, 1, sharex=True, layout='tight')
    ax = np.atleast_1d(ax)

    if isinstance(title, str):
        ax[0].set_title(title)

    if semilogx:
        ax[0].set_xscale('log')

    plot_curves(ax[0], f, h_mag, linestyles, colors)

    ax[0].set_ylabel(f'Magnitude [{scale}]')
    ax[0].grid(which='both')
    ax[0].format_coord = lambda x, y: f'x={x:.3f}, y={y:.3f}'

    if phase:
        
        plot_curves(ax[1], f, np.angle(h), linestyles, colors)
        
        ax[1].set_ylabel('Phase [radians]')
        ax[1].set_xlabel('Frequency [Hz]')
        ax[1].set_yticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
        ax[1].set_yticklabels(['$-\pi$', '$-\pi/2$', '0', '$\pi/2$', '$\pi$'])
        ax[1].set_xlim([f[0], f[-1]])
        ax[1].format_coord = lambda x, y: f'x={x:.3f}, y={y:.2f}'
        ax[1].grid(which='both')
    
    ax[-1].set_xticks([1, 20, 50, 100, 200, 500, 1e3, 2e3, 5e3, 10e3, 20e3, 100e3, 1e6])
    ax[-1].set_xticklabels(['1', '20', '50', '100', '200', '500', '1k', '2k', '5k', '10k', '20k', '100k', '1M'])
    ax[-1].set_xlim([f[0], f[-1]])

    plot_legend(ax, h, legend, colors, linestyles)
    plt.show()


def plot_curves(ax, x, y, linestyles, colors):
    if y.ndim == 2:
        for j, y_array in enumerate(y):
            ax.plot(x, y_array, color=colors[j % len(colors)])
    else:
        for i, y_slice in enumerate(y):
            for j, y_array in enumerate(y_slice):
                ax.plot(x, y_array, linestyles[i % len(linestyles)], color=colors[j % len(colors)])
    
    

def plot_legend(ax, h, legend, colors, linestyles):

    colors_nbr, linestyles_nbr = len(colors), len(linestyles)

    if legend is not None: 

        if isinstance(legend, str):
            ax[0].plot([], [], '-', label=legend)
            
        else:
            if h.ndim == 3:
                first_key = list(legend.keys())[0]
                first_key_values = legend[first_key]

                for i in range(len(first_key_values)):
                    ax[0].plot([], [], linestyles[i % linestyles_nbr], color='grey', label=f'{first_key}: {first_key_values[i]}')
    
            last_key = list(legend.keys())[-1]  #if h.ndim==2, last_key is obviously the only key in the legend dictionary
            last_key_values = legend[last_key]

            for i in range(len(last_key_values)):
                ax[0].plot([], [], '-', color=colors[i % colors_nbr], label=f'{last_key}: {last_key_values[i]}')
        ax[0].legend(loc='best')


