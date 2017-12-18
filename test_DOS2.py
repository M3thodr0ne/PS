execfile("DOS_alter.py")

execfile("Files.py")

os.chdir("Results16alter111")

T = readIteration(str(9))
well_energies = T[0]
z = T[3]
Ef = T[1]
V = T[4]

a = 3.9*10**(-10)
size = 20
param = 0.0
ntot = 0.25

acc = 0.0001



#Ef = fermiEnergyComputer(well_energies,param,size,a,acc,ntot)

print(Ef)

#Ef = 0.2957

size2 = 10

Lkx = [[] for i in range(4*len(well_energies)*len(well_energies[0]))]
Lky = [[] for i in range(4*len(well_energies)*len(well_energies[0]))]

for i in range(2*size2):
	for j in range(2*size2):
		k = Kmesh(a,2*size2,size2-i,size2-j)
		B = band_struct(well_energies,param,k,a)
		for index in range(len(B)):
			if abs(Ef-real(B[index])) < 0.01:
				Lkx[index].append(k[0]*a/pi)
				Lky[index].append(k[1]*a/pi)

#for i in range(len(Lkx)):
	#plot(Lkx[i],Lky[i],'rx')

#plt.show()

Pdens = bandPopComputer(well_energies,param,size,a,Ef)

P = rho3d_alter(z,Ef,well_energies,V,Pdens)

plot(z,P)
plt.show()
