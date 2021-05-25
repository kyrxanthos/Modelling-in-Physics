"""
Modelling and Visualisations in Physics
Mar/2021
Checkpoint 3a: Cahn Hilliard PDE Solver
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
    parser.add_argument("--sweeps", type=int, help="sweeps for the simulation", default=100000)
    parser.add_argument("--dim", type=int, help="Dimensions of lattice", default=100)
    parser.add_argument("--kappa", type=float, help="positive constant related to surface tension", default=0.1)
    parser.add_argument("--a", type=float, help="constant ensuring stability", default=0.1)
    parser.add_argument("--M", type=float, help="mixing constant", default=0.1)
    parser.add_argument("--delta_x", type=float, help="spatial discretisation", default=1.5)
    parser.add_argument("--delta_t", type=float, help="time discretisation", default=2.0)
    parser.add_argument("--phi0", type=float, help="order parameter", default=0.)
    parser.add_argument("--measurement", type=str, help="kind of measurement", default="animation", choices=["animation", "free_energy"])

    args = parser.parse_args()
    return args


class CahnHilliard():

    def __init__(self, sweeps, dim, phi0, kappa, a, delta_x, delta_t, M, measurement):
        self.sweeps = sweeps
        self.dim = dim
        self.phi_latt = np.empty((self.dim, self.dim))
        self.mu_latt = np.empty((self.dim, self.dim))
        self.phi0 = phi0
        self.get_phi()
        self.kappa = kappa
        self.a = a
        self.delta_x = delta_x
        self.delta_t = delta_t
        self.M = M
        self.measurement = measurement

    def get_phi(self):
        # asssigns random noise and intial phi0 to phi.
        for i in range(self.dim):
            for j in range(self.dim):
                self.phi_latt[i, j] = self.phi0 + np.random.uniform(-0.1, 0.1)
    

    def update(self):
        #update function with discretisation for the CH Equation
        kappa = self.kappa
        a = self.a
        delta_x = self.delta_x
        delta_t = self.delta_t
        M = self.M

        mu_new = -a * self.phi_latt + a * self.phi_latt ** 3 - \
            (kappa / (delta_x)**2) * (np.roll(self.phi_latt, -1, axis=0) + np.roll(self.phi_latt, 1, axis=0) - 2 * self.phi_latt \
                + np.roll(self.phi_latt, -1, axis=1) + np.roll(self.phi_latt, 1, axis=1) -2 * self.phi_latt)

        phi_new = self.phi_latt + ((M * delta_t) / (delta_x) ** 2) * (np.roll(mu_new, -1, axis=0) + \
            np.roll(mu_new, 1, axis=0)  + np.roll(mu_new, -1, axis=1) + np.roll(mu_new, 1, axis=1) - 4 * mu_new )

        self.phi_latt = phi_new

    def Grad(self):
        # returns x,y components of the gradient
        grad = np.gradient(self.phi_latt)
        return np.array(grad)[0], np.array(grad)[1]


    def free_energy(self):
        # caclulates the free energy of the system.
        FE = 0.
        a = self.a
        kappa = self.kappa

        FE = - (a / 2) * self.phi_latt ** 2 + \
            (a / 4) * self.phi_latt ** 4 + \
            (kappa / (2 * (self.delta_x) ** 2)) * (self.Grad()[0] ** 2 + self.Grad()[1] ** 2)

        return np.sum(FE)


    def sim(self):
        #simulation of evolution of CH Equation
        times = []
        FEs = []
        for n in range(self.sweeps):
            # times.append(n)
            if n % 100 == 0 :
                if self.measurement == 'animation':
                    plt.cla()
                    im = plt.imshow(self.phi_latt, vmin=-1., vmax=1., cmap='PiYG', animated=True)
                    plt.pause(0.0001)

                print('sweep {}'.format(n))
                times.append(n)

                FE = self.free_energy()
                FEs.append(FE)

            self.update()

        return times, FEs

    def plots(self,x,y, outdir):
        plt.figure(figsize = (10,8))
        plt.title(r'Free Energy Density against Sweeps')
        plt.xlabel(r"Sweeps")
        plt.ylabel(r"Free Energy Density")
        plt.plot(x,y)
        plt.savefig(f'{outdir}/FE_plot_dim_{self.dim}_phi_{self.phi0}.png')

        
def main():

    args = parse_args()

    sweeps= args.sweeps
    dim = args.dim
    kappa = args.kappa
    a = args.a
    M = args.M 
    delta_x = args.delta_x
    delta_t = args.delta_t
    phi0 = args.phi0
    measurement = args.measurement

    #create output folder
    if not os.path.isdir("Results"):
        os.makedirs("Results")
    outdir = 'Results'
    
    CH = CahnHilliard(sweeps, dim, phi0, kappa, a, delta_x, delta_t, M, measurement)
    times, FEs = CH.sim()

    with open(f"{outdir}/FE_data_{phi0}.dat", 'w+') as f:
        for i, j in zip(times, FEs):
            f.write("{0:f} {1:f} \n".format(i, j))

    CH.plots(times, FEs, outdir)

main()
