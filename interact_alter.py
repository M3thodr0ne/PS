execfile("DOS_alter_interact.py")


def interaction_loop_alter(well_energies,ntot,Ef,U,Nstep,param,size,a,acc):

    BD = bandPopComputer(well_energies,param,size,a,Ef)
    print(sum(BD))
    Efcorr = Ef

    for i in range(Nstep):
        Efmin0 = Efcorr - 0.25*U
        Efmax0 = Efcorr + 0.25*U
        Efcorr = fermiEnergyComputerInter(well_energies,param,size,a,acc,ntot,Efmin0,Efmax0,U,BD)

        BD2 = bandPopComputerInter(well_energies,param,size,a,Efcorr,U,BD)
        BD = [[0.6*BD[i][j] + 0.4*BD2[i][j] for j in range(4*len(well_energies[i]))] for i in range(len(well_energies))]
        for l in range(3):
            print(sum(BD[l]))

    return [BD,Efcorr]
