import sympy as sp
from component_class import *



class Netlist:    
    def __init__(self, inputList):

        self.components         = [self.create_component(item) for item in inputList]
        self.passiveComponents  = self.get_components_of_type(self.components, PassiveComponent)
        self.voltageSources     = self.get_components_of_type(self.components, VoltageSource)
        self.currentSources     = self.get_components_of_type(self.components, CurrentSource)

        self.nodeNbr = self.get_nodes_nbr()

        self.A = self.create_A_matrix()
        self.x = self.create_x_vector()
        self.b = self.create_b_vector()


    def create_component(self, netlist_line):
        start_node, end_node, symbol, value = netlist_line.split()

        if symbol.startswith('V'):
            return VoltageSource(start_node, end_node, symbol, value)
        
        elif symbol.startswith('I'):
            return CurrentSource(start_node, end_node, symbol, value)
        
        elif symbol.startswith('R'):
            return Resistance(start_node, end_node, symbol, value)
        
        elif symbol.startswith('C'):
            return Capacitor(start_node, end_node, symbol, value)
        
        elif symbol.startswith('L'):
            return Inductance(start_node, end_node, symbol, value)
        
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
    


    def create_admittance_matrix(self):
        n = self.nodeNbr
        Y = sp.zeros(n, n)

        for component in self.passiveComponents:
            # Convert node names to indexed integer for matrix
            
            i = int(component.start_node)
            j = int(component.end_node)

            admittance = component.admittance

            # Diagonal elements 
            #Each element in the diagonal matrix is equal to the sum of the admittance 
            #of each element connected to the corresponding node. So the first diagonal element is the sum 
            #of admittance connected to node 1, the second diagonal element is the sum of conductances 
            #connected to node 2, and so on.

            Y[i, i] += admittance
            Y[j, j] += admittance
            
            # Off-diagonal elements

            Y[i, j] -= admittance
            Y[j, i] -= admittance

        #return the "real" admittance matrix, since the first row and first column are 
        #not included in the real one (the node 0 which refer to the ground is not considered)
        return Y[1:,1:]


    def create_B_matrix(self):
        
        n = self.nodeNbr 
        m = len(self.voltageSources)

        B = sp.zeros(n , m)
        
        i = 0
        
        for source in self.voltageSources:
            B[int(source.start_node), i] = -1
            B[int(source.end_node),   i] = 1
            i+=1
        
        return B[1:,::]
        

    def create_A_matrix(self):
        Y = self.create_admittance_matrix()
        B = self.create_B_matrix()
        
        A = sp.BlockMatrix([[Y, B], 
                            [B.T, sp.zeros(len(B[0,:]), len(B[0,:]))]])
        
        return A.as_explicit()


    def create_x_vector(self):
        n = self.nodeNbr - 1

        #iV_vector holds the unknown currents through the voltage sources
        #v_vector holds the unknown voltages at each node, except at the ground

        v_vector = sp.Matrix([sp.symbols(f'v{j+1}') for j in range(n)])
        iV_vector = sp.Matrix([sp.symbols(f'i{source.symbol}') for source in self.voltageSources])
        x_vector = v_vector.col_join(iV_vector)
        
        return x_vector
    

    def create_b_vector(self):
        
        #z_vector is structured as a column vector of the form [i_vector, e_vector]
        #The i_vector is a vector for which each index corresponding to a particular node. 
        #The value of each element of i_vector is determined by the sum of current sources into the corresponding node. 
        #If there are no current sources connected to the node, the value is zero. The node 0 is not considered.
        #e_vector holds the known voltage of each voltage sources

        i_vector = sp.zeros(self.nodeNbr, 1)

        for source in self.currentSources:
            i_vector[int(source.start_node)] -= source.symbol
            i_vector[int(source.end_node)] += source.symbol

        
        i_vector =  i_vector[1:,:] #avoid the node 0

        e_vector = sp.Matrix([source.symbol for source in self.voltageSources])
        
        b_vector = i_vector.col_join(e_vector)
        
        return b_vector
    

    def display_components(self, components_list = None):

        if components_list is None:
            components_list = self.components

        print('-'*73)
        print(f" {'Type':<14} | {'Symbol':^6} | {'Start':^5} | {'End':^5} | {'Admittance':^15} | {'Value':^11} ")
        print('-'*73)
        for component in components_list:
            print(f" {str(component):<14} | {str(component.symbol):^6} | {component.start_node:^5} | {component.end_node:^5} | {str(component.admittance):^15} | {component.value:^11} ")
        print('-'*73)




    
    


    
    
    
    
    
    
    
    
