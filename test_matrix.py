execfile("DOS_alter.py")

well_energies = [[i for i in range(2)] for j in range(3)]

a = 3.9*10**(-10)
size = 40
param = 0.03
ntot = 0.24336

M = matrixBands(well_energies,param,[0,0],a)

for i in range(len(M)):
    print(M[i])
