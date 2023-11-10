# -*- coding: utf-8 -*-
"""
Created on Sun Oct 22 14:25:17 2023

@author: eliot
"""
from component_class import VoltageSource, CurrentSource, Capacitor, Resistance, Inductance
import sympy as sp

def create_component(netlist_line):
    start_node, end_node, symbol = netlist_line.split()

    if symbol.startswith('V'):
        return VoltageSource(start_node, end_node, symbol)
    
    elif symbol.startswith('I'):
        return CurrentSource(start_node, end_node, symbol)
    
    elif symbol.startswith('R'):
        return Resistance(start_node, end_node, symbol)
    
    elif symbol.startswith('C'):
        return Capacitor(start_node, end_node, symbol)
    
    elif symbol.startswith('L'):
        return Inductance(start_node, end_node, symbol)
     
    else:
        raise ValueError(f"Unknown component symbol: {symbol}")


def display_components(components_list):
    for component in components_list:
        print(f" {component} {component.symbol}: start node {component.start_node},end node {component.end_node}, admittance {component.admittance}")

    
def get_admittance(component):
    """
    Returns the admittance of the component.
    """
    
    if isinstance(component, Resistance) or isinstance(component, Capacitor) or isinstance(component, Inductance):
        return component.admittance
    else:
        return 0

def get_nodes_nbr(components):
    """
    Returns the total number of unique nodes from the list of components. 
    The ground node is included 
    The set data type is used for avoid to count multiple time the same node (no duplication)
    """
    nodes = set()
    for component in components:
        nodes.add(component.start_node)
        nodes.add(component.end_node)
    return len(nodes)


def get_components_of_type(components, component_type):
    """
    Returns the total number of one type of component from the list of components. 
    The arg component_type must be a subclass of the Component Class
    """
    components_list = [component for component in components if isinstance(component, component_type)]
    return components_list
        


def create_admittance_matrix(components):
    n = get_nodes_nbr(components)
    Y = sp.zeros(n, n)

    for component in components:
        # Convert node names to indexed integer for matrix
        
        i = int(component.start_node)
        j = int(component.end_node)
        admittance = get_admittance(component)


        # Diagonal elements 
        #Each element in the diagonal matrix is equal to the sum of the admittance 
        #of each element connected to the corresponding node. So the first diagonal element is the sum 
        #of admittance connected to node 1, the second diagonal element is the sum of conductances 
        #connected to node 2, and so on.
        Y[i, i] += admittance
        Y[j, j] += admittance
        
        # Off-diagonal elements
        if i != j:
            Y[i, j] -= admittance
            Y[j, i] -= admittance

    #return the "real" admittance matrix, since the first row and first column are 
    #not included in the real one (the node 0 which refer to the ground is not considered)
    return Y[1:,1:]


def create_B_matrix(components):
    sources_list = get_components_of_type(components, VoltageSource)
    
    n = get_nodes_nbr(components)
    m = len(sources_list)

    B = sp.zeros(n , m)
    
    i = 0
    
    for source in sources_list:
        B[int(source.start_node), i] = -1
        B[int(source.end_node),   i] = 1
        i+=1
    
    return B[1:,::]
    

def create_A_matrix(components):
    Y = create_admittance_matrix(components)
    B = create_B_matrix(components)
    
    A = sp.BlockMatrix([[Y, B], 
                        [B.T, sp.zeros(len(B[0,:]), len(B[0,:]))]])
    
    return A.as_explicit()


def create_x_vector(components):
    sources_list = get_components_of_type(components, VoltageSource)
    n = get_nodes_nbr(components)-1

    #iV_vector holds the unknown currents through the voltage sources
    #v_vector holds the unknown voltages at each node, except at the ground

    iV_vector = sp.Matrix([sp.symbols(f'i{source.symbol}') for source in sources_list])
    v_vector = sp.Matrix([sp.symbols(f'v{j+1}') for j in range(n)])
    x_vector = v_vector.col_join(iV_vector)
    
    return x_vector
    

def create_z_vector(components):
    
    #z_vector is structured as a column vector if the form [i_vector, e_vector]
    
    
    #The i_vector is a vector for which each index corresponding to a particular node. 
    #The value of each element of i_vector is determined by the sum of current sources into the corresponding node. 
    #If there are no current sources connected to the node, the value is zero.
    #The node 0 is not considered.
    
    
    #e_vector holds the known voltage of each voltage sources

    c_sources_list = get_components_of_type(components, CurrentSource)
    v_sources_list = get_components_of_type(components, VoltageSource)
    
    n = get_nodes_nbr(components)
    
    i_vector = sp.zeros(n, 1)

    i = 0
    
    for source in c_sources_list:
        i_vector[int(source.start_node)] += -source.symbol
        i_vector[int(source.end_node)] += source.symbol
        i+=1

    
    i_vector =  i_vector[1:,:] #avoid the node 0

    e_vector = sp.Matrix([source.symbol for source in v_sources_list])
    
    z_vector = i_vector.col_join(e_vector)
    
    return z_vector
                 

def solve_x_vector(netlist):
   #return the voltage at each node of the circuit.
   #and current passing through each voltage generator
    
    components_list = [create_component(item) for item in netlist]
    
    A = create_A_matrix(components_list)
    x = create_x_vector(components_list)
    z = create_z_vector(components_list)
    
    solution = A.LUsolve(z)
    
    return sp.simplify(solution)




    
    


    
    
    
    
    
    
    
    