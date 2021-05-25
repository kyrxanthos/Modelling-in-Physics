"""
Modelling and Visualisations in Physics
Feb/2021
Checkpoint 2a: Simulation of Game of Life
Kyriacos Xanthos

"""


import os
import math
import random
import matplotlib.animation as animation
import numpy as np
import argparse
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
    parser.add_argument("--init", type=str, help="Initial conditions of model", default="random", choices=["random", "oscilator", "glider", "absorbing"])
    parser.add_argument("--sweeps", type=int, help="sweeps for the simulation", default=10000)
    parser.add_argument("--dim", type=int, help="Dimensions of lattice", default=50)
    parser.add_argument("--measurement", type=str, help="Kind of measurement", default="animation", choices=["animation", "equilibration", "glider_speed"])
    args = parser.parse_args()
    return args


class Game_of_Life():

    def __init__(self, lx, ly, config = 'random'):
        self.lx = lx
        self.ly = ly
        self.init = config
        self.eq = False
        if config == 'random':
            self.latt = np.random.choice([1, 0], size=(self.lx, self.ly))

        elif config  ==  'absorbing':
            self.latt = np.zeros((self.lx, self.ly), dtype=float)
            self.latt[int(self.lx / 2 - 1),int(self.ly /2)] = 1
            self.latt[int(self.lx / 2 + 2),int(self.ly /2)] = 1
            self.latt[int(self.lx / 2),int(self.ly /2 + 1)] = 1
            self.latt[int(self.lx / 2),int(self.ly /2 - 1)] = 1
            self.latt[int(self.lx / 2 + 1),int(self.ly /2 + 1)] = 1
            self.latt[int(self.lx / 2 + 1),int(self.ly /2 - 1)] = 1

        elif config == 'oscilator':
            self.latt = np.zeros((self.lx, self.ly), dtype=float)
            self.latt[int(self.lx / 2),int(self.ly /2)] = 1
            self.latt[int(self.lx / 2 ), int(self.ly / 2 + 1)] = 1
            self.latt[int(self.lx / 2 ), int(self.ly / 2 - 1)] = 1

        elif config == 'glider':
            self.latt = np.zeros((self.lx, self.ly), dtype=float)
            self.latt[int(self.lx / 2 - 1),int( self.ly /2)] = 1
            self.latt[int(self.lx / 2),int( self.ly /2 + 1)] = 1
            self.latt[int(self.lx / 2 + 1),int( self.ly /2 - 1)] = 1
            self.latt[int(self.lx / 2 + 1),int( self.ly /2 + 1)] = 1
            self.latt[int(self.lx / 2 + 1),int( self.ly /2)] = 1



    def nearest_neighbours(self, cell):
        # returns the number of alive neighbours
        alive_cells = 0
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

        # 8 nearest neighbours
        nearest_neighbours_list = [self.latt[iup, j], self.latt[iup, jup], self.latt[iup, jdown],
                              self.latt[idown, j], self.latt[idown, jup], self.latt[idown, jdown],
                              self.latt[i,  jup], self.latt[i,jdown]]

        for cell in nearest_neighbours_list:
            if cell == 1:
                alive_cells += 1

        return alive_cells


    def evolve(self):
        # evolves the simulation with the rules for the dynamics
        updated_latt = np.copy(self.latt)
        for i in range(self.lx):
            for j in range(self.ly):
                alive = self.nearest_neighbours([i,j])
                if updated_latt[i,j] == 1:
                    if alive < 2 or alive > 3:
                        updated_latt[i,j] = 0
                elif updated_latt[i,j] == 0 and alive == 3:
                    updated_latt[i, j] = 1
        return updated_latt

    def sim(self, latt, sweeps):
        # runs the simulatioon
        for _ in range(sweeps):
            plt.cla()
            plt.title(f'Game of life {self.lx}x{self.ly} (condition: {self.init})')
            plt.imshow(latt.latt, animated = True)
            plt.draw()
            plt.pause(0.0001)            
            latt.latt =latt.evolve()


    def get_live(self):
        #counts live cells
        return np.sum(self.latt)

    def steady_state(self, cells):
        # Boolean returning True if the number of total alive cells 
        # is the same for the last 4 times we ran the simulation
        if len(cells) > 4:
            if cells[-1] == cells[-2] and cells[-2] == cells[-3] and cells[-3] == cells[-4]:
                self.eq =  True


    def check_boundary(self):
        # checks if the glider is not crossing any boundaries
        alive_cell = np.where(self.latt == 1)
        xpos = alive_cell[0] ; ypos = alive_cell[1]
        if np.amax(xpos) - np.amin(xpos) <= 2 and np.amax(ypos) - np.amin(ypos) <= 2:
            return True
        else:
            return False

    def glider_position(self):
        # gets the x and y positions of the glider
        alive_cell = np.where(self.latt == 1)
        xpos = alive_cell[0] ; ypos = alive_cell[1] 
        return xpos, ypos

    def glider_com(self, xs, ys):
        # calculates the centre of mass given the x and y positions of the glider
        x_com = np.sum(xs) * (1 / (self.get_live()))
        y_com = np.sum(ys) * (1 / (self.get_live()))
        return x_com, y_com

    def plot_hist(self, times_list, outdir):
        #plots histogram
        plt.figure(figsize=(10,6))
        plt.title("Time needed for simulations to reach Equilibrium")
        plt.xlabel("Time (Sweeps)")
        plt.ylabel("Frequency")
        n, bins, _ = plt.hist(times_list, bins = np.arange(0, 3000, 100))
        plt.savefig(f'{outdir}/histogram_dim_{self.lx}.png')
        return n, bins

    def plot_glider_speed(self, times, xpos, ypos, zpos, cut, outdir):
        # fits the x,y,r using polyfit and plots the results
        p1 = np.polyfit(times[:cut], xpos[:cut], 1)
        y_fit1 = np.array(times[:cut])*p1[0] + p1[1]

        p2 = np.polyfit(times[:cut], ypos[:cut], 1)
        y_fit2 = np.array(times[:cut])*p2[0] + p2[1]

        p3 = np.polyfit(times[:cut], zpos[:cut], 1)
        y_fit3 = np.array(times[:cut])*p3[0] + p3[1]

        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize = (18,6))
        fig.suptitle(f'Glider speed {self.lx}x{self.ly}')
        ax1.plot(times[:cut], y_fit1, label = f"x(t) = {p1[0]:.3}t + {p1[1]:.3}")
        ax1.plot(times[:cut], xpos[:cut], "rx")
        ax1.set_xlabel('time')
        ax1.set_ylabel('x(t)')
        ax1.legend()

        ax2.plot(times[:cut], y_fit2, label = f"y(t) = {p2[0]:.3}t + {p2[1]:.3}")
        ax2.plot(times[:cut], ypos[:cut], "rx")
        ax2.set_xlabel('time')
        ax2.set_ylabel('y(t)')
        ax2.legend()

        ax3.plot(times[:cut], y_fit3, label = f"r(t) = {p3[0]:.3}t + {p3[1]:.3}")
        ax3.plot(times[:cut], zpos[:cut], "rx")
        ax3.set_xlabel('time')
        ax3.set_ylabel('r(t)')
        ax3.legend()

        plt.tight_layout()

        fig.savefig(f'{outdir}/glider_speed_plots_dim_{self.lx}.png')



def main():
    args = parse_args()

    # Directory to store all the results
    if not os.path.isdir("Results"):
        os.makedirs("Results")
    outdir = 'Results'

    init = args.init
    sweeps = args.sweeps
    dim = args.dim
    measurement = args.measurement

    if measurement == 'animation':
        # performs animation
        GOL = Game_of_Life(dim, dim, init)
        GOL.sim(GOL, sweeps)

    if measurement == 'equilibration':
        # equilibrium point at which this stops evolving, 
        # you should plot a distribution of the times needed to equilibrate.
        f = open(f"{outdir}/{measurement}_data_dim_{dim}.dat","w+")

        'sim to check active sites'
        times = []
        for sweep in range(sweeps):
            print(f"sweep {sweep} out of {sweeps}")
            alive_cells = []
            GOL = Game_of_Life(dim,dim, 'random')
            n_sweeps = 0
            while GOL.eq == False:
                GOL.latt = GOL.evolve()
                n_sweeps += 1
                # gets all live neighbours 
                alive_tot = GOL.get_live()
                alive_cells.append(alive_tot)
                GOL.steady_state(alive_cells)
                # so that we do not get in an infinite loop.
                if n_sweeps == 5000:
                    break
            # subtract the 4 that have reached equilibration
            times.append(len(alive_cells) - 4 )
            f.write("{0:f}\n".format(len(alive_cells) - 4 ))


        n, bins = GOL.plot_hist(times, outdir)

        with open(f"{outdir}/histogram_data_dim_{dim}.dat", 'w+') as f:
            for i, j in zip(n, bins):
                f.write("{0:f} {1:f} \n".format(j, i))


    if measurement == 'glider_speed':

        GOL = Game_of_Life(dim,dim, 'glider')
        xpos =  [] ; ypos =  [] ; times = []

        for sweep in range(sweeps):
            if sweep % 100 == 0:
                print(f"sweep {sweep} out of {sweeps}")

            GOL.latt = GOL.evolve()
            # if glider is not at the boundary, record position and centre of mass
            if  GOL.check_boundary():
                xs , ys = GOL.glider_position()
                xcom, ycom = GOL.glider_com(xs, ys)
                xpos.append(xcom) ; ypos.append (ycom)
                times.append(sweep)
        zpos = np.sqrt(np.square(xpos) + np.square(ypos))

        GOL.plot_glider_speed(times, xpos, ypos, zpos, 15, outdir)

        with open(f"{outdir}/{measurement}_data.dat", 'w+') as f:
            for i, j, k, l in zip(times, xpos, ypos, zpos):
                f.write("{0:f} {1:f} {2:f} {3:f} \n".format(i, j, k, l))


main()





    
    