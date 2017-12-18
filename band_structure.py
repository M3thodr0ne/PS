import os

execfile("DOS_alter_interact.py")
execfile("Files.py")



def hexagPath(resol,a):
    pathx = [0 for i in range(2*resol)]
    pathy = [0 for i in range(2*resol)]
    labelnames = ["K","$\Gamma$","M"]
    for i in range(resol):
        #K-G path
        k = Kmesh(a,resol,(resol-i),(resol-i))
        pathx[i] = 0.5*(resol-i)
        pathy[i] = 0.*(resol-i)
        #G-M path
        k = Kmesh(a,resol,i,0)
        pathx[i+resol] = 0.33*i
        pathy[i+resol] = -0.33*i
    return [pathx,pathy]

def plotBands(well_energ,param,a,resol,U,Mpopbands,Ef):
    P = hexagPath(resol,a)
    pathx = P[0]
    pathy = P[1]
    tabEn = [[0 for j in range(len(pathx))] for i in range(len(well_energ)*len(well_energ[0])*4)]
    for i in range(len(pathx)):
        k = Kmesh(a,resol,pathx[i],pathy[i])
        energ = band_struct_Inter(well_energ,param,k,a,U,Mpopbands)
        #energ = band_struct(well_energ,param,k,a)
        energ = sort(real(energ))
        for j in range(len(energ)):
            tabEn[j][i] = real(energ[j]) - Ef + 0. * j
    for j in range(len(energ)):
        plot(tabEn[j])
    plot([0 for i in range(len(pathx))],"--",color="black")
    ylim(-0.2,0.1)
    xlim(resol/2,3*resol/2)
    plt.show()
    return tabEn

def transfMatrix(well_energ,tabEn,tabColor):
    Nbands = len(well_energ[0])
    tabinterm = [[(tabEn[i][j],tabColor[i][j]) for i in range(len(tabEn))] for j in range(len(tabEn[0]))]
    tabEnsort = [[0 for j in range(len(tabEn[0]))] for i in range(len(tabEn))]
    for i in range(len(tabinterm)):
        M = sorted(tabinterm[i],key=lambda colonnes: colonnes[1])
        for j in range(len(M)):
            tabEnsort[j][i] = M[j][0]
    return tabEnsort


def plotBandsimproved(well_energ,param,a,resol,U,Mpopbands,Ef):
    P = hexagPath(resol,a)
    pathx = P[0]
    pathy = P[1]
    Nbands = len(well_energ[0])
    tabEn = [[0 for j in range(len(pathx))] for i in range(3*Nbands*4)]
    tabColor = [[0 for j in range(len(pathx))] for i in range(3*Nbands*4)]
    for i in range(len(pathx)):
        k = Kmesh(a,resol,pathx[i],pathy[i])
        B = eig(matrixBandsInteract(well_energ,param,k,a,U,Mpopbands))
        energies = B[0]
        eigenvect = real(B[1])
        for j in range(len(energies)):
            k = 0
            colo = "black"
            while k < Nbands:
                if sum(abs(eigenvect[ 4*Nbands*l + k * 4 + m ][j])**2 for l in range(3) for m in range(4)) > 10**(-1):
                    indexband = k
                    tabColor[j][i] = k
                k += 1
            tabEn[j][i] = real(energies[j]) - Ef + 0. * j
    tSort = transfMatrix(well_energ,tabEn,tabColor)
    tSort2 = [[0 for j in range(len(pathx))] for i in range(3*Nbands*4)]
    for i in range(len(tSort[0])):
        for j in range(Nbands):
            tband = sort([tSort[12*j + m][i] for m in range(12)])
            for w in range(len(tband)):
                tSort2[w+j*12][i] = tband[w]
    for j in range(Nbands):
        for w in range(12):
            if w == 0:
                plot([tSort2[w+12*j][i] for i in range(len(tSort2[0]))],color = colorBand(j),label = str(j))
            else:
                plot([tSort2[w+12*j][i] for i in range(len(tSort2[0]))],color = colorBand(j))
    plot([0 for i in range(len(pathx))],"--",color="black")
    ylim(-0.2,0.1)
    xlim(resol/2,3*resol/2)
    legend()
    plt.show()

def fileBand(filename):
    TT = readIteration(str(filename))
    well_energies = TT[0]
    Ef = TT[1]
    a = 3.9*10**(-10)
    resol = 15
    param = 0.03
    U = 2.7*0.
    Mpopbands = bandPopComputer(well_energies,param,resol,a,Ef)
    plotBands(well_energies,param,a,resol,U,Mpopbands,Ef)

def colorBand(i):
    if i == 0:
        return "blue"
    if i == 1:
        return "red"
    if i == 2:
        return "yellow"
    if i == 3:
        return "green"
    else:
        return "black"


chemin = os.getcwd()
makeFold("111",chemin)
chemin = os.getcwd()
makeFold("Kmesh15Rashba",chemin)


TT = readIteration(str(0))
well_energies = TT[0]
Ef = TT[1]
a = 3.9*10**(-10)
resol = 15
param = 0.03
U = 2.7*1.
acc = 0.0001
ntot = 0.25
Mpopbands = bandPopComputer(well_energies,param,resol,a,Ef)

#fileBand(12)

print(Ef)

Eftild = fermiEnergyComputerInter(well_energies,param,resol,a,acc,ntot,Ef-0.2,Ef+ntot*U+0.2,U,Mpopbands)

print(Eftild)

plotBandsimproved(well_energies,param,a,resol,U,Mpopbands,Eftild)
