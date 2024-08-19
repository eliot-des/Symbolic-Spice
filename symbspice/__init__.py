import sympy as sp
import numpy as np
from .submodules.TransferFunction import TransferFunction
from .submodules.components import AdmittanceComponent, Resistance, Capacitor, Inductance, VoltageSource, ExternalVoltageSource, CurrentSource, IdealOPA, Transformer, Gyrator 

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
        self.n = self.__get_nodes()   # Number of nodes including the ground node

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
    
    def __get_nodes(self):
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
            component.stamp_MNA()

        #set the symbolic value at each voltage node in the x 'state' vector:
        for i in range(self.n):
            self.x[i] = sp.symbols(f'v{i}')

        #slice the A matrix, b and x vectors to remove the datum node (ground node)
        self.A = self.A[1:, 1:]
        self.b = self.b[1:, :]
        self.x = self.x[1:, :]

    def show(self):
        print('\n\nA matrix:\n')
        sp.pprint(self.A)
        print('\n\nx vector:\n')
        sp.pprint(self.x)
        print('\n\nb vector:\n')
        sp.pprint(self.b)

    def __solve(self, simplify = False, use_symengine = False):
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


    def tf(self, output_node, input_node=1):
        """
        Returns transfer function object of the netlist.

        Parameters
        ----------
        output_node (int): Node where the output is set.
        input_node (int): Node where Vin is set.
        norm (bool): Normalize coefficients (b/a0 & a/a0 if True).
        """
        
        
        # Ensure the system is solved:
        if not hasattr(self, 'x_solution'):
            self.__solve()

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
            node = np.argwhere(netlist[:,0] == 'Vin')[0]
        elif np.any(netlist[:, 0] == 'V1'):
            node = np.argwhere(netlist[:,0] == 'V1')[0]
        else:
            raise Exception("Netlist is missing Vin, add a voltage source named Vin where the input is wanted")
        
        # Make sure it's value is 1
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
    
        # If there is a net label for the output set it as the last node
        if np.any(netlist[:, 1:3] == 'out'):
            # Replace out to max node
            netlist[:,1:3] = np.char.replace(netlist[:,1:3], 'out', str(self.outNode))

        # return netlist as list of rows (each row is a component: [name, start node, end node, value])
        outlist = []
        for n in range(netlist.shape[0]):
            # Clean non-numeric values such as R, C, etc
            if np.char.isdigit(netlist[n,3][0]) == False:
                netlist[n,3] = '0'
                
            outlist.append(' '.join(str(item) for item in netlist[n,:]))

        return outlist