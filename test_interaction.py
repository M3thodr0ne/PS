import os

execfile("Schrodinger.py")
execfile("Poisson.py")
execfile("DOS.py")
execfile("Files.py")



def weightEigen(name):
    T = readIteration(name)
    well_energies = T[0]
    Ef = T[1]
    n2d = T[2]
    z = T[3]
    V = T[4]
    tabweight = [0 for i in range(len(well_energies[0]))]

    for i in range(len(well_energies[0])):
        psi = ith_eigenfuncacc(well_energies[0][i],mz[0],V,z)
        plot(z,[abs(psi[k])**4 for k in range(len(psi))])
        tabweight[i] = sum(abs(psi[k])**4 for k in range(len(psi)))
    plt.show()
    return tabweight


chemin = os.getcwd()
makeFold("111",chemin)
chemin = os.getcwd()
makeFold("Kmesh15Rashba",chemin)

TT = readIteration(str(0))
well_energies = TT[0]
Ef = TT[1]
a = 3.9*10**(-10)
resol = 30
param = 0
resolE = 20

Tab = weightEigen("0")

print(Tab)
