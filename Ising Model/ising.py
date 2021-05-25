"""
Modelling and Visualisations in Physics
Feb/2021
Checkpoint 1: Simulation of Ising Model
Kyriacos Xanthos

"""

import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['grid.alpha'] = 0.5
plt.rc('grid', linestyle="--", color='grey')
matplotlib.pyplot.title(r'ABC123 vs $\mathrm{ABC123}^{123}$')
import os
import sys
import math
import random
import matplotlib.animation as animation
import numpy as np
import argparse


def parse_args():
    parser = argparse.ArgumentParser() 
    parser.add_argument("--dynamics", type=str, help="Dynamics of Ising model", default="Glauber", choices=["Glauber", "Kawasaki"])
    parser.add_argument("--sweeps", type=int, help="sweeps for the simulation", default=10000)
    parser.add_argument("--dim", type=int, help="Dimensions of lattice", default=50)
    parser.add_argument("--animate", type=bool, help="Show Animation", default=False)
    parser.add_argument("--temp", type=float, help="Single temperature for the animation", default=1.2, required=False)
    args = parser.parse_args()
    return args


class ising_spin():

    def __init__(self, lx, ly, config = 'random'):
        self.lx = lx
        self.ly = ly
        if config == 'random':
            self.spin = 2*np.random.randint(2, size=(lx,ly))-1
        elif config == 'up':
            self.spin = np.ones((self.lx, self.ly), dtype=float)
        elif config == 'half':
            x = np.zeros((self.lx,int(self.lx/2))) -1
            y = np.ones((self.ly,int(self.ly/2)))
            self.spin = np.concatenate((y,x), axis=1)

    
    def E_calc(self):
        # Calculates the total energy of the whole lattice
        E = 0. 
        for i in range(self.lx):
            for j in range(self.ly):
                iup = i + 1
                if (i == self.lx - 1): 
                    iup = 0
                jup = j + 1
                if (j == self.ly - 1):
                    jup = 0
                E += - self.spin[i,j] * (self.spin[iup, j] + self.spin[i, jup])
        return E

    def E_calc_nearest(self, sp):
        # Energy Calculation for nearest neighbours
        E = 0.
        iup = sp[0] + 1 ; jup = sp[1] + 1
        idown = sp[0] - 1 ; jdown = sp[1] - 1
        # Periodic boundary conditions
        if (sp[0] == self.lx -1):
            iup = 0
        if (sp[0] == 0):
            idown = self.lx -1
        if (sp[1] == self.ly -1):
            jup = 0
        if (sp[1] == 0):
            jdown = 0
        E += - self.spin[sp[0],sp[1]] * (self.spin[iup, sp[1]] + self.spin[sp[0], jup] + \
            self.spin[idown, sp[1]] + self.spin[sp[0], jdown])
        return E


    def random_int(self):
        # picks a random position for a spin
        x = random.randint(0, self.lx - 1)
        y = random.randint(0, self.ly - 1)
        return [x,y]

    def check_adjacency(self, sp1, sp2):
        # checks if two spins are adjacent
        iup = sp1[0] + 1 ; jup = sp1[1] + 1 
        idown = sp1[0] - 1 ; jdown = sp1[1] - 1 
        # B.C.
        if (sp1[0] == self.lx - 1):
            iup = 0
        if (sp1[1] == self.ly -1):
            idown = 0
        if (sp1[0] == 0):
            idown = self.lx - 1
        if (sp1[1] == 0):
            jdown = self.ly - 1
        if iup == sp2[0] or idown == sp2[0] or jup == sp2[1] or jdown == sp2[1] :
            return True
        else:
            return False
        
    def switch(self, sp1, sp2):
        # Switch the spins with each other
        dummy_1 = self.spin[sp1[0], sp1[1]]
        dummy_2 = self.spin[sp2[0], sp2[1]]
        self.spin[sp2[0], sp2[1]] = dummy_1
        self.spin[sp1[0], sp1[1]] = dummy_2


    def glauber_metropolis(self, T):
        # Metropolis Test for glauber
        flip_attempt = self.random_int()
        Delta_E = -2 * self.E_calc_nearest(flip_attempt) 
        if Delta_E <= 0. :
            #flip spin
            self.spin[flip_attempt[0], flip_attempt[1]] *= -1
        elif Delta_E > 0. :
            prob = math.exp(- Delta_E / T)
            p = np.minimum(1, prob)
            if random.random() < p :
                #flip spin
                self.spin[flip_attempt[0], flip_attempt[1]] *= -1

    def kawasaki_metropolis(self, T):
        # Metroopolis Test for Kawasaki
        Delta_E = 0 
        flip_attempt1 = self.random_int()
        flip_attempt2 = self.random_int()
        # check if the random spins are the same
        while(self.spin[flip_attempt1[0], flip_attempt1[1]] == self.spin[flip_attempt2[0], flip_attempt2[1]]):
            flip_attempt1 = self.random_int()
            flip_attempt2 = self.random_int()
        # calculation for nearest neighbours
        if self.check_adjacency(flip_attempt1, flip_attempt2):
            Delta_E = -2 * self.E_calc_nearest(flip_attempt1) \
                -2 * self.E_calc_nearest(flip_attempt2) + 4. 
        # calculation for not nearest neighbours
        else:
            Delta_E = -2 * self.E_calc_nearest(flip_attempt1) \
                -2 * self.E_calc_nearest(flip_attempt2)
        # now choose to make the switch or not
        if Delta_E <= 0.:
            #always switch
            self.switch(flip_attempt1, flip_attempt2)
        else:
            prob = math.exp(- Delta_E / T)
            p = np.minimum(1, prob)
            if random.random() < p :
                #switch spins
                self.switch(flip_attempt1, flip_attempt2)

    def magnetisation(self):
        # Calculates the magnetisation
        M = 0.
        for i in range(self.lx):
            for j in range(self.ly):
                M += self.spin[i, j]
        return abs(M)

    def normal_error(self, lst):
        # Calculates standard error of an array
        err = np.sqrt((np.mean(np.square(lst)) \
            - np.mean(lst) ** 2)/(len(lst) - 1))
        return err


    def bootstrap(self, lst, T, k, kind):
        # Calculates error using the Bootstrap method for either energies or magnetisations
        chi_average = chi_averagesq = 0.
        C_average = C_averagesq = 0.
        for i in range(k):
            valtot = valtotsq = 0.
            for element in range(len(lst)):
                index = random.randint(0, len(lst) - 1)
                val = lst[index]
                valtot += val
                valtotsq += val * val
            if kind == 'magnetisation':
                chi = (1 / (self.lx * self.ly * T)) * \
                    ((valtotsq/len(lst)) - (valtot/len(lst))**2)
                chi_average += chi
                chi_averagesq += chi * chi
            elif kind =='energy':
                C = (1 / (self.lx * self.ly * T**2)) * \
                    ((valtotsq/len(lst)) - (valtot/len(lst))**2)
                C_average += C
                C_averagesq += C * C

        if kind == 'magnetisation': 
            error = math.sqrt(abs((chi_averagesq / k) - (abs(chi_average / k)) ** 2))
        if kind == 'energy':
            error = math.sqrt(abs((C_averagesq / k) - (abs(C_average / k)) ** 2))
        return error



    def run(self, spin, sweeps, T, animate, dynamics, threshold):
        # Run dynamics of ising model
        E_all =[] ; M_all = []
        count = 0
        # <E>, <E^2>, <M>, <M^2>
        Et = Mt = Etsq = Mtsq = 0. 
        for sweep in range(sweeps):
            for point in range(self.lx * self.ly):
                if dynamics == 'Glauber':
                    spin.glauber_metropolis(T)
                elif dynamics == 'Kawasaki':
                    spin.kawasaki_metropolis(T)
            # Take uncorrelated measurements
            if sweep % 10 == 0:
                # Wait for equilibration
                if sweep >= threshold:
                    # Keep count to divide total measurements and find mean
                    count += 1
                    E =spin.E_calc()
                    M = spin.magnetisation()
                    E_all.append(E) ; M_all.append(M)
                    Et += E ; Mt += M
                    Etsq += E**2 ; Mtsq += M**2

            if animate == True:
                plt.cla()
                im = plt.imshow(spin.spin, vmin=-1., vmax=1., cmap='jet', animated=True)
                plt.draw()
                plt.pause(0.0001)

        
        chi = (1 / (self.lx * self.ly * T)) * ((abs(Mtsq/count)) - (abs(Mt/count))**2)
        C = (1 / (self.lx * self.ly * T**2)) * ((Etsq/count) - (Et/count)**2)
        E = Et / count
        M = Mt / count

        return E, M, C, chi, spin, E_all, M_all

    def plots(self, E, M, chi, C, T, dynamics, sweeps, E_error, M_error, chi_error, C_error, outdir):
        # Plots of energy, magnetisation, susceptibility and specific heat against temperature
        fig = plt.figure(figsize=(18, 10))

        ax =  fig.add_subplot(2, 2, 1 )
        plt.errorbar(T, E, yerr=E_error, fmt = '.', color ='red')
        plt.xlabel("Temperature (T)", fontsize=20)
        plt.ylabel("Energy ", fontsize=20);         plt.axis('tight')

        ax =  fig.add_subplot(2, 2, 2 )
        plt.errorbar(T, M, yerr=M_error, fmt = '.', color ='blue')
        plt.xlabel("Temperature (T)", fontsize=20) 
        plt.ylabel("Magnetization ", fontsize=20);   plt.axis('tight')

        ax =  fig.add_subplot(2, 2, 3 )
        plt.errorbar(T, C, yerr=C_error, fmt = '.', color ='red')
        plt.xlabel("Temperature (T)", fontsize=20)  
        plt.ylabel("Specific Heat ", fontsize=20);   plt.axis('tight')  

        ax =  fig.add_subplot(2, 2, 4 )
        plt.errorbar(T, chi, yerr=chi_error, fmt = '.', color ='blue')
        plt.xlabel("Temperature (T)", fontsize=20) 
        plt.ylabel("Susceptibility", fontsize=20);   plt.axis('tight')
        fig.suptitle(r'{} Dynamics, lattice dimensions: {}$\times${}'.format(dynamics, self.lx, self.ly), fontsize = 30)
        plt.savefig("{}/plots_{}_{}.png".format(outdir, dynamics, sweeps), dpi=fig.dpi)



def main():
    args = parse_args()

    # Directory to store all the results
    if not os.path.isdir("Results"):
        os.makedirs("Results")
    outdir = 'Results'

    T = args.temp
    lx = args.dim
    ly = lx
    sweeps = args.sweeps
    dynamics = args.dynamics
    # number of temperature points
    nt = 21
    # array of temperatures to check
    temperatures = np.linspace(1,3, nt) 
    threshold = 100

    # Initiate lists to store results
    E_mean_list, M_mean_list, chi_list, C_list = np.zeros(nt), np.zeros(nt), np.zeros(nt), np.zeros(nt)
    E_error_list, M_error_list, chi_error_list, C_error_list  = np.zeros(nt), np.zeros(nt), np.zeros(nt), np.zeros(nt)

    # Select initial configuration depending on Dynamics
    if args.dynamics == 'Glauber':
        config = 'up'
    elif args.dynamics == 'Kawasaki':
        config = 'half'
    else:
        config = 'random'

    # Show animation, chose random inital state to show equilibration
    if args.animate:
        spin = ising_spin(lx, ly, config = 'random')
        spin.run(spin, sweeps, T, args.animate, dynamics, threshold)

    # Calculation for all temperature points
    if not args.animate:
        spin = ising_spin(lx, ly, config)
        f = open("{}/measurements_{}_{}_{}.dat".format(outdir, args.dynamics, args.sweeps, args.dim),"w+")

        for n, T in enumerate(temperatures):
            # Track progression
            if n % 2 == 0:
                print('Temperature point {} out of {}'.format(n, len(temperatures)))
            # Run simulation
            E_mean, M_mean, C, chi, spin, E_all, M_all = spin.run(spin, sweeps, T, args.animate, dynamics, threshold)

            # Store measurements
            E_mean_list[n] = E_mean
            M_mean_list[n] = np.abs(M_mean)
            chi_list[n] = chi
            C_list[n] = C

            # Calculate and store errors
            E_error = spin.normal_error(E_all)
            M_error = spin.normal_error(M_all)
            chi_error = spin.bootstrap(M_all, T, k = 1000, kind = 'magnetisation')
            C_error = spin.bootstrap(E_all,T,k=1000, kind = 'energy')
            E_error_list[n] = E_error
            M_error_list[n] = M_error
            chi_error_list[n] = chi_error
            C_error_list[n] = C_error

            # Output measurements in a .dat file for each temperature
            f.write('{0:f} {1:f} {2:f} {3:f} {4:f} {5:f} {6:f}\n'.format(T, E_mean, C , M_mean, chi, chi_error, C_error))

        # Print critical temperature as predicted
        ind1 = np.unravel_index(np.argmax(chi_list, axis=None), chi_list.shape)
        print('Temperature at the maximum of susceptibility: {}'.format(temperatures[ind1]))
        ind2 = np.unravel_index(np.argmax(C_list, axis=None), C_list.shape)
        print('Temperature at the maximum of specific heat: {}'.format(temperatures[ind2]))

        # Plot results
        spin.plots(E_mean_list, M_mean_list, chi_list, C_list, 
        temperatures, dynamics, sweeps, E_error_list, M_error_list, 
        chi_error_list, C_error_list, outdir)


main()





    
    