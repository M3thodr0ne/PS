execfile("DOS_alter.py")

execfile("Files.py")

well_energies = [[i for i in range(1)] for j in range(3)]

a = 3.9*10**(-10)
size = 40
param = 0.0
ntot = 0.24336

acc = 0.00001



Ef = fermiEnergyComputer(well_energies,param,size,a,acc,ntot)


print(Ef)

#Ef = 0.2957

size2 = 40

Lkx = [[] for i in range(4*len(well_energies)*len(well_energies[0]))]
Lky = [[] for i in range(4*len(well_energies)*len(well_energies[0]))]

for i in range(2*size2):
	for j in range(2*size2):
		k = Kmesh(a,2*size2,(size2-i),(size2-j))
		B = band_struct(well_energies,param,k,a)
		for index in range(len(B)):
			if abs(Ef-real(B[index])) < 0.005:
				Lkx[index].append(k[0]*a/pi)
				Lky[index].append(k[1]*a/pi)

for i in range(len(Lkx)):
	plot(Lkx[i],Lky[i],'rx')

B = eig(matrixBands(well_energies,param,k,a))
eigenval = B[0]
eigenvect = B[1]

print(real(eigenval))

plt.show()

Pdens = bandPopComputer(well_energies,param,size,a,Ef)

print(Pdens)
print(sum(Pdens))
