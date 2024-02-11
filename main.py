from netlist_class import *

#declaration of the netlist
#arbitrary values
inputList  =   ['Vin 0 1 1', 
                'R1 1 2 8',
                'L1 2 3 4',
                'C1 3 0 1']


netlist = Netlist(inputList)

netlist.display_components()


#system of the form A·x = b
A = netlist.A
x = netlist.x
z = netlist.b

solution =  sp.simplify(A.LUsolve(z))


print('\n\nA matrix :\n ')
sp.pprint(A)
print('\n\nx vector :\n ')
sp.pprint(x)
print('\n\nb vector :\n ')
sp.pprint(b)
print('\n\nSolution :\n ')
sp.pprint(sp.Eq(x , solution))

