import os

os.chdir("/Users/Pierre/Desktop/Stage_4A/Poisson_Schrodinger/Simulations/src/")

execfile("veigen_plot.py")

number_file = 19

os.chdir("/Users/Pierre/Desktop/Stage_4A/Poisson_Schrodinger/Simulations/Results10/")

plotDensityMulti(str(number_file))

os.chdir("/Users/Pierre/Desktop/Stage_4A/Poisson_Schrodinger/Simulations/Results15/")

plotDensityMulti(str(number_file))


os.chdir("/Users/Pierre/Desktop/Stage_4A/Poisson_Schrodinger/Simulations/Results16/")

plotDensityMulti(str(number_file))


os.chdir("/Users/Pierre/Desktop/Stage_4A/Poisson_Schrodinger/Simulations/Results17/")

plotDensityMulti(str(number_file))

os.chdir("/Users/Pierre/Desktop/Stage_4A/Poisson_Schrodinger/Simulations/Results20/")

plotDensityMulti(str(number_file))


plt.legend()
plt.show()