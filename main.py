from netlist_class import *

#declaration of the netlist
inputList  =   ['0 1 Vin 1', 
                '1 2 R1 8',
                '2 3 L1 4',
                '2 0 R2 1',
                '2 0 L2 1',
                '3 0 C1 1']


netlist = Netlist(inputList)

netlist.display_components()


#system of the form Ax = b
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

