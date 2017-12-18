import os

execfile("maps_alter.py")



TT = readIteration(str(0))
well_energies = TT[0]
Ef = TT[1]
a = 3.9*10**(-10)
resol = 20
param = 0
resolE = 20

#compute the relative weight below some energy and make a color plot

T = mapOrbitalComputer(well_energies,param,resol,a,Ef)
mkx = T[0]
mky = T[1]
pxy = T[2]
pyz = T[3]
pxz = T[4]

x = np.linspace(-1,1,2*resol)
y = np.linspace(-1,1,2*resol)

X , Y = np.meshgrid(x,y)

z = np.array([pxy[i][j] for i in range(2*resol) for j in range(2*resol)])


Z = z.reshape(2*resol,2*resol)



plt.pcolor(X,Y,Z)

plt.colorbar()
plt.show()

z = np.array([pxz[i][j] for i in range(2*resol) for j in range(2*resol)])


Z = z.reshape(2*resol,2*resol)



plt.pcolor(X,Y,Z)

plt.colorbar()
plt.show()

z = np.array([pyz[i][j] for i in range(2*resol) for j in range(2*resol)])


Z = z.reshape(2*resol,2*resol)



plt.pcolor(X,Y,Z)

plt.colorbar()
plt.show()

#make a plot of the DOS

T2 = mapDOS(well_energies,param,resol,resolE,a,Ef)

TE = T2[0]
DE = T2[1]
plot(TE,DE)

plt.show()
