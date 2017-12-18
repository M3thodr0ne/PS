execfile("DOS_alter_interact.py")

execfile("Files.py")

def mapOrbitalComputer(well_energies,param,size,a,Ef):
    pxy = [[0 for i in range(2*size)] for j in range(2*size)]
    pxz = [[0 for i in range(2*size)] for j in range(2*size)]
    pyz = [[0 for i in range(2*size)] for j in range(2*size)]
    meshpointx = [[0 for i in range(2*size)] for j in range(2*size)]
    meshpointy = [[0 for i in range(2*size)] for j in range(2*size)]
    for i in range(2*size):
        for j in range(2*size):
            k = Kmesh(a,2*size,(-size+i),(-size+j))
            meshpointx[i][j] = k[0]
            meshpointy[i][j] = k[1]
            B = eig(matrixBands(well_energies,param,k,a))
            eigenval = B[0]
            eigenvect = B[1]
            N = len(eigenval)/3
            for index in range(len(eigenval)):
                if eigenval[index] < Ef:
                    pxy[i][j] = sum(abs(eigenvect[      index2][index])**2 for index2 in range(N))
                    pxz[i][j] = sum(abs(eigenvect[N   + index2][index])**2 for index2 in range(N))
                    pyz[i][j] = sum(abs(eigenvect[2*N + index2][index])**2 for index2 in range(N))

    return [meshpointx,meshpointy,pxy,pyz,pxz]

def mapDOS(well_energies,param,sizek,size,a,Ef):
    TE = [Ef-0.5 - Ef/size * i for i in range(size)]
    DOSE = [0 for i in range(size)]
    for i in range(2*sizek):
        for j in range(2*sizek):
            k = Kmesh(a,2*size,(-size+i),(-size+j))
            B = eig(matrixBands(well_energies,param,k,a))
            eigenval = B[0]
            for index in range(len(eigenval)):
                indexE = int(-size/Ef * (real(eigenval[index])-Ef+0.5))
                if indexE < len(TE):
                    DOSE[indexE] += 1
    return [TE,DOSE]

def mapDOSorb(well_energies,param,sizek,size,a,Ef):
    TE = [Ef-0.5 + 1./size * i for i in range(size)]
    DOSxy = [0 for i in range(size)]
    DOSxz = [0 for i in range(size)]
    DOSyz = [1 for i in range(size)]
    for i in range(2*sizek):
        for j in range(2*sizek):
            k = Kmesh(a,2*size,(-size+i),(-size+j))
            B = eig(matrixBands(well_energies,param,k,a))
            eigenval = B[0]
            eigenvect = B[1]
            N = len(eigenval)/3
            for index in range(len(eigenval)):
                indexE = int(-size/Ef * (real(eigenval[index])-Ef+1.0))
                if indexE < len(TE):
                    DOSxy[indexE] += sum(abs(eigenvect[      index2,index])**2 for index2 in range(N))
                    DOSxz[indexE] += sum(abs(eigenvect[N   + index2,index])**2 for index2 in range(N))
                    DOSyz[indexE] += sum(abs(eigenvect[2*N + index2,index])**2 for index2 in range(N))

    return [TE,DOSxy,DOSxz,DOSyz]


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


T2 = mapDOSorb(well_energies,param,resol,resolE,a,Ef)

TE = T2[0]
DOSxy = T2[1]
DOSxz = T2[2]
DOSyz = T2[3]

plot(TE,DOSxy,label="xy")
plot(TE,DOSxz,label="xz")
plot(TE,DOSyz,label="yz")
plt.axvline(Ef,color = "black",linestyle = '--')

plt.legend()
plt.show()
