import os

execfile("interact_alter.py")
execfile("Files.py")



chemin = os.getcwd()
makeFold("111",chemin)
chemin = os.getcwd()
makeFold("Kmesh15test",chemin)

U = 3
U1 = 0

T = readIteration(str(0))
well_energies = T[0]
Ef = T[1]
ntotal = T[2]

Nstep = 30
param = 0
size = 15

a = 3.9*10**(-10)
acc = 0.0001
k = [0,0]

Mpopbands = bandPopComputer2(well_energies,param,size,a,Ef)
#
#print(matrixBandsInteract(well_energies,param,k,a,U,Mpopbands))
#print(matrixBandsInteract(well_energies,param,k,a,0,Mpopbands))

print(Mpopbands)
print(sum(Mpopbands))
print(ntotal)

#T = interaction_loop_alter(well_energies,ntotal,Ef,U,Nstep,param,size,a,acc)

#for i in range(len(T)):
	#print(sum(T[i]))
