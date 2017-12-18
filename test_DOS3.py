execfile("DOS_alter.py")

execfile("Files.py")



os.chdir("Results16alter111")

T = readIteration(str(9))
well_energies = T[0]
z = T[3]
Ef = T[1]
V = T[4]

a = 3.9*10**(-10)
size = 15
param = 0.0
ntot = 0.25

Tabacc = [10**(-(2+0.5*i)) for i in range(6)]

tabEf = [0 for i in range(len(Tabacc))]
tabN = [0 for i in range(len(Tabacc))]

for i in range(len(Tabacc)):
    tabEf[i] = fermiEnergyComputer(well_energies,param,size,a,Tabacc[i],ntot)
    tabN[i] = compute_density(well_energies,param,tabEf[i],size,a)

semilogx(Tabacc,tabN)
plt.title("Size of K-mesh = "+str(size)+"x"+str(size))
plt.xlabel("Accuracy of Fermi energy determination")
plt.ylabel("Density")
namepic = "Accuracy_size"+str(size)+".pdf"
savePicture(namepic)

plt.show()

semilogx(Tabacc,tabEf)
plt.title("Size of K-mesh = "+str(size)+"x"+str(size))
plt.xlabel("Accuracy of Fermi energy determination")
plt.ylabel("Fermi Energy")
namepic = "Accuracy_size_Ef"+str(size)+".pdf"

savePicture(namepic)

plt.show()
