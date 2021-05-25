"""
Modelling and Visualisations in Physics
Feb/2021
Checkpoint 2b: Simulation of SIRS
Kyriacos Xanthos

"""

import os
import math
import random
import numpy as np
import argparse
import time
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size': 13})
plt.rcParams['grid.alpha'] = 0.5
plt.rc('grid', linestyle="--", color='grey')
matplotlib.pyplot.title(r'ABC123 vs $\mathrm{ABC123}^{123}$')

def parse_args():
    parser = argparse.ArgumentParser() 
    parser.add_argument("--sweeps", type=int, help="sweeps for the simulation", default=1000)
    parser.add_argument("--dim", type=int, help="Dimensions of lattice", default=50)
    parser.add_argument("--p1", type=float, help="Probability of S -> I", required= False)
    parser.add_argument("--p2", type=float, help="Probability of I -> R", required= False)
    parser.add_argument("--p3", type=float, help="Probability of R -> S", required= False)
    parser.add_argument("--measurement", type=str, help="Kind of measurement", default="animation", \
        choices=["animation", "phase_space", "variance", "vaccine"])

    args = parser.parse_args()
    return args


class SIRS():

    def __init__(self, lx, ly, p1, p2, p3, config = 'random'):
        self.lx = lx
        self.ly = ly
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        if config == 'random':
            self.latt = np.random.choice([0, 1, 2], size=(self.lx, self.ly))
        

    def nearest_neighbours(self, cell):
        # checks if at least one nearest neighbour is infected
        i =  cell[0] ; j = cell[1]
        iup = i + 1 ; jup = j + 1
        idown = i - 1 ; jdown = j - 1
        # Periodic boundary conditions
        if (i == self.lx -1):
            iup = 0
        if (i == 0):
            idown = self.lx -1
        if (j == self.ly -1):
            jup = 0
        if (j == 0):
            jdown = self.ly -1

        nearest_neighbours_list = [self.latt[iup, j], self.latt[idown, j],self.latt[i,  jup], self.latt[i,jdown]]

        for cell in nearest_neighbours_list:
            if cell == 1:
                return True
        return False

    def evolve(self):
        # evolves sirs simulation with the set of rules
        for _ in range(self.lx * self.ly):
            # picks random index (monte carlo) 
            index = np.random.choice(self.latt.shape[0], 2, replace=True)  
            cell = self.latt[index[0], index[1]]
            rand = random.random()
            if cell == 0 and self.nearest_neighbours(index) == True:
                if rand < self.p1:
                    cell = 1
            elif cell == 1:
                if rand < self.p2:
                    cell = 2
            elif cell == 2:
                if rand < self.p3:
                    cell = 0
            else:
                cell = cell
            # updates the cell
            self.latt[index[0], index[1]] = cell

    def get_infected(self):
        # measures the total infected in the lattice
        tot_inf = 0.
        for i in range(self.lx):
            for j in range(self.ly):
                if self.latt[i,j] == 1:
                    tot_inf += 1
        tot_inf
        return tot_inf

    def sim(self, system, sweeps):
        # runs the simulation
        for _ in range(sweeps):
            plt.cla()
            plt.title(f'SIRS Simulation for p1 = {self.p1}, p2 = {self.p2}, p3 = {self.p3}')
            plt.imshow(system.latt, vmin=0., vmax=2., cmap='jet', animated=True)
            plt.draw()
            plt.pause(0.0001)
            system.evolve()


    def bootstrap(self, lst, k):
        # Calculates error using the Bootstrap method
        # lst - list to perform resampling
        # k - number of resamplings 
        observable_average = observable_averagesq = 0.
        for _ in range(k):
            valtot = valtotsq = 0.
            for _ in range(len(lst)):
                index = random.randint(0, len(lst) - 1)
                val = lst[index]
                valtot += val
                valtotsq += val * val

            # observable = (valtotsq - valtot**2 ) / (self.lx * self.ly)
            observable = ((valtotsq) / (len(lst)) - ((valtot) / (len(lst)))**2) / (self.lx * self.ly)
            observable_average += observable
            observable_averagesq += observable * observable

        error = math.sqrt(abs((observable_averagesq / k) - (abs(observable_average / k)) ** 2))

        return error

    def phase_plotter(self, space, outdir, kind):
        # plots phase space for average infected mean and variance
        plt.figure(figsize = (10,10))
        plt.title(f'Average Infected ({kind})')
        plt.xlabel(r"$p3$")
        plt.ylabel(r"$p1$")
        im = plt.imshow(space, cmap="hot", origin='lower', extent=(0, 1, 0, 1))
        plt.colorbar()
        plt.savefig("{}/phase_diagram_{}.png".format(outdir, kind))
        plt.clf()  

    def plot_variance(self, x, y, error, outdir):
        # plots variance against p1
        plt.figure(figsize = (10,8))
        plt.title(f'Variance of the number of infected sites, p3 = {self.p3}')
        plt.xlabel(r"$p1$")
        plt.ylabel(r"$\frac{<{I}^2> - {<I>}^2}{N}$")
        plt.errorbar(x, y, yerr = error, fmt='.b')
        plt.savefig(f'{outdir}/variance_plot_dim{self.lx}_p3_{self.p3}.png')

    def plot_vaccine(self, x, y, error, outdir):
        # plots infected fraction with immunity
        plt.figure(figsize = (10,8))
        plt.title(r'Average infected fraction versus fraction of imunity $f_{Im}$')
        plt.xlabel(r"$f_{Im}$")
        plt.ylabel(r"$\frac{<I>}{N}$")
        plt.errorbar(x, y, yerr = error, fmt='.b')
        plt.savefig(f'{outdir}/vaccine_plot_dim{self.lx}_p3_{self.p3}.png')
                    

def main():
    args = parse_args()

    sweeps = args.sweeps
    dim = args.dim
    measurement = args.measurement
    # threshold for equilibration
    threshold = 100

    p1 = args.p1
    p2 = args.p2
    p3 = args.p3
    lx = ly = dim

    if not os.path.isdir("Results"):
        os.makedirs("Results")
    outdir = 'Results'

    t1 = time.time()


    if measurement == 'animation':
        system = SIRS(dim, dim, p1, p2, p3, config= 'random')
        system.sim(system, sweeps)

    if measurement == 'phase_space':

        p1s = np.linspace(0,1, 21)
        p3s = np.linspace(0,1, 21)
        p2 = 0.5

        phase_space = np.zeros((len(p1s), len(p3s)))
        variance_space = np.zeros((len(p1s), len(p3s)))

        for i, p1 in enumerate(p1s):
            print('p1 {} out of {}'.format(i, len(p1s)))
            # store data per p1
            inf_p1 = []
            inf_p1_var = []
            for p3 in p3s:
                # store data per p3
                inf_p3 = []
                system = SIRS(lx, ly, p1, p2, p3 , config= 'random')
                for sweep in range(sweeps):
                    system.evolve()
                    #equilibration
                    if sweep >= threshold:
                        current_infected = system.get_infected()
                        if current_infected == 0.:
                            # terminatewhen simulation reaches an absorbing state
                            break
                        else:
                            inf_p3.append(current_infected)
                
                # ensure list has been populated before taking mean/var
                if len(inf_p3) != 0:
                    inf_p1.append(np.mean(inf_p3) / (dim * dim))
                    inf_p1_var.append(np.var(inf_p3) / (dim * dim))

                else:
                    inf_p1.append(0.)
                    inf_p1_var.append(0.)

            phase_space[:,i ] = inf_p1
            variance_space[:, i] = inf_p1_var


        with open(f"{outdir}/{measurement}_data.dat", 'w+') as f:
            for i, j, k, l in zip(p1s, p3s, inf_p1, inf_p1_var):
                f.write("{0:f} {1:f} {2:f} {3:f} \n".format(i, j, k, l))

        system.phase_plotter(phase_space, outdir,  kind = 'Mean')
        system.phase_plotter(variance_space, outdir,  kind = 'Variance')

    if measurement == 'variance':
        p1s = np.linspace(0.2,0.5, 51)
        # fixed p3 and p2
        p3 = 0.5
        p2 = 0.5

        var_array = np.zeros(len(p1s))
        error_array = np.zeros(len(p1s))

        for i, p1 in enumerate(p1s):
            print(f"{i} out of {len(p1s)}")
            tot_inf = []
            system = SIRS(lx, ly, p1, p2, p3 , config= 'random')
            for sweep in range(sweeps):
                system.evolve()
                if sweep >= threshold:
                    tot_inf.append(system.get_infected())

            # variance
            # var_array[i] = np.var(tot_inf) / (lx * ly)
            t_i = np.asarray(tot_inf)
            var_array[i] = (np.mean(t_i**2) - (np.mean(t_i))**2) / (lx * ly)
            # error using resampling method
            error_array[i] = system.bootstrap(tot_inf, 100)


        with open(f"{outdir}/{measurement}_data.dat", 'w+') as f:
            for i, j, k in zip(p1s, var_array, error_array):
                f.write("{0:f} {1:f} {2:f} \n".format(i, j, k))
            

        system.plot_variance(p1s, var_array, error_array, outdir)

    if measurement == 'vaccine':
        fIm  = np.linspace(0, 0.5 ,51)

        all_infected = []
        errors = []

        for i in range(5):
            # list of infected for each i
            inf_i = []
            print(f"{i} out of 5")

            for f in fIm:
                # print(f)
                # list of infected for each f
                inf_f = []
                system = SIRS(lx, ly, p1, p2, p3 , config= 'random')

                # fraction of immune f for lattice lx * ly
                for _ in range(int(f * lx * ly)):
                    # choose a random cell
                    index = np.random.choice(system.latt.shape[0], 2, replace=True)  
                    # change that cell to a state '3' which stays unchanged 
                    # throughout simulation (immune).
                    system.latt[index[0], index[1]] = 3

                for sweep in range(sweeps):
                    system.evolve()
                    #equilibration
                    if sweep >= threshold:
                        if system.get_infected() == 0.:
                            # total infected 0 means disease is eradicated so
                            # we can stop simulation
                            if len(inf_f) == 0:
                                inf_f.append(0.)
                            inf_f[:] = [0.] * len(inf_f)
                            break
                        else:
                            inf_f.append(system.get_infected() / (lx * ly))

                inf_i.append(np.mean(inf_f))

            all_infected.append(inf_i)

        inf_av = np.mean(all_infected, axis = 0)
        # reshape to have a list for calculating the error
        inf_reshaped = np.asarray(all_infected).T

        for j in inf_reshaped:
            err = np.std(j) / np.sqrt(len(j))
            errors.append(err)

        with open(f"{outdir}/{measurement}_dim{lx}_data.dat", 'w+') as f:
            for i, j, k in zip(fIm, inf_av, errors):
                f.write("{0:f} {1:f} {2:f} \n".format(i, j, k))

        system.plot_vaccine(fIm, inf_av, errors, outdir)

    t2 = time.time()
    
    print (f"Solution calculated in {t2-t1} seconds.")

main()

