execfile("DOS_alter.py")

execfile("Files.py")



os.chdir("Results16alter111")

T = readIteration(str(9))
well_energies = T[0]
z = T[3]
Ef = T[1]
V = T[4]

a = 3.9*10**(-10)
size = 10
param = 0.0
ntot = 0.25

acc = 10**(-5)

T = fermiEnergyComputerTot(well_energies,param,size,a,acc,ntot)

Tabacc = [T[0][-i-1] for i in range(10)]
tabEf = [T[1][-i-1] for i in range(10)]
tabN = [T[2][-i-1] for i in range(10)]


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
