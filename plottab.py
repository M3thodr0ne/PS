#this file to plot the evolution of some quantities with respect to the density of the 2DEG

import os
import numpy as np
from pylab import *
from scipy import *
from math import *
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import fsolve



#density of the 2DEG
tabn2d = [1.0, 1.5, 1.55, 1.6, 1.7, 2.0]
tabdensity = [tt*10**18*(3.9*10**(-10))**2 for tt in tabn2d]

#mean electric field felt by the 2DEG
tabEmoy = [63, 208, 229, 251, 297, 470]

#Fermi energy
tabEf = [254, 269, 285, 315, 417]

#energy of the dxy1
tabdxy1 = [-112, -118, -124, -135, -173]

#energy of the dxy2
tabdxy2 = [-51, -53, -55, -59, -72]

#energy of the dxz1
tabdxz1 = [-31,-33,-33,-35,-41]

#energy of the crossing
tabCross = [-10.5,-10.5,-10,-9.2,-6.2]

plot(tabdensity,tabEmoy)

xlabel("Density (electron per UC)",fontsize = 16)
ylabel("Mean electric field (MV.m$^{-1}$)", fontsize = 16)

legend(loc = 3)
plt.show()