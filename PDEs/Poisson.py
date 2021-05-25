"""
Modelling and Visualisations in Physics
Mar/2021
Checkpoint 3b: Poisson PDE Solver 3D
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
    parser.add_argument("--dim", type=int, help="Dimensions of lattice", default=50)
    parser.add_argument("--phi0", type=float, help="order parameter", default=0.)
    parser.add_argument("--tolerance", type=float, help="tolerance to determine convergence", default=0.01)
    parser.add_argument("--omega", type=float, help="relaxation parameter", default=1.9)
    parser.add_argument("--field_type", type=str, help="Magnetic or Electric Field", default="E", choices=["E", "M"])
    parser.add_argument("--SOR_calc", type=bool, help="Find optimum omega", default=False)
    parser.add_argument("--algorithm", type=str, help="Algorithm to solve Poissons Equation", default="jacobi", \
        choices=["jacobi", "gs", "gs_SOR"])

    args = parser.parse_args()
    return args

class poisson():

    def __init__(self, dim, phi0):
        self.dim = dim
        self.phi0 = phi0
        self.phi_latt = np.empty((self.dim, self.dim, self.dim))
        self.get_phi()
        # monopole
        self.rho = np.zeros((self.dim, self.dim, self.dim))
        self.rho[int(self.dim / 2)][int(self.dim / 2)][int(self.dim / 2)] = 1.


    def get_phi(self):
        # initiate phi with appropriate B.C.
        for i in range(self.dim):
            for j in range(self.dim):
                for k in range(self.dim):
                    if i == 0 or j == 0 or k == 0 or i == self.dim - 1 or j == self.dim - 1 or k == self.dim - 1:
                        # Dirichlet B.C.
                        self.phi_latt[i,j,k] = 0.
                    else:
                        # update phi with noise.
                        # self.phi_latt[i, j, k] = self.phi0 + np.random.uniform(-0.1, 0.1)
                        self.phi_latt[i, j, k] = self.phi0 


    def nearest_neighbours(self, latt, i, j, k):
        # checks if at least one nearest neighbour is infected
        iup = i + 1 ; jup = j + 1 ; kup = k + 1
        idown = i - 1 ; jdown = j - 1 ; kdown = k - 1

        # Periodic boundary conditions
        if (i == self.dim -1):
            iup = 0
        if (i == 0):
            idown = self.dim -1
        if (j == self.dim -1):
            jup = 0
        if (j == 0):
            jdown = self.dim -1
        if (k == self.dim - 1):
            kup = 0
        if (k == 0):
            kdown = self.dim - 1

        nearest_neighbours_list = [latt[iup, j, k], latt[idown, j, k], latt[i,  jup, k], \
            latt[i,jdown, k], latt[i, j, kup], latt[i, j, kdown]]

        return nearest_neighbours_list

    def boundaries(self, img, width=1):
        #sets boundaries to zero
        img[:width, :, :] = 0
        img[-width:, :, :] = 0
        img[:, :width, :] = 0
        img[:, -width:, :] = 0
        img[:, :, -width:] = 0 
        img[: , :, :width] = 0
        return img


    def jacobi(self):
        phi_new = (1/6) * (self.boundaries(np.roll(self.phi_latt, -1, axis=0)) + 
        self.boundaries(np.roll(self.phi_latt, 1, axis=0)) +
        self.boundaries(np.roll(self.phi_latt, -1, axis=1)) +
        self.boundaries(np.roll(self.phi_latt, 1, axis=1)) +
        self.boundaries(np.roll(self.phi_latt, -1, axis=2)) +
        self.boundaries(np.roll(self.phi_latt, 1, axis=2)) +
        self.rho)

        return phi_new

    def gauss_seidel(self, omega, SOR):
        #gs
        phi_new = np.copy(self.phi_latt)
        for i in range(self.dim):
            for j in range(self.dim):
                for k in range(self.dim):
                    nn = self.nearest_neighbours(phi_new, i, j, k)

                    if i == 0 or j == 0 or k == 0:
                        phi_new[i, j, k] = 0.
                    else:
                        if SOR == True:
                            phi_new[i, j, k] = (1/6) * (sum(nn) + self.rho[i,j,k]) * omega + (1 - omega) * phi_new[i,j,k]
                        else:
                            phi_new[i, j, k] = (1/6) * (sum(nn) + self.rho[i,j,k])
        return phi_new


    def Efield(self):
        # returns x,y,z components of E field
        E = np.gradient(self.phi_latt)
        return -1*np.array(E)[0], -1*np.array(E)[1], -1*np.array(E)[2]

    def Mfield(self):
        # curl
        M = np.gradient(self.phi_latt)
        Mx = M[1] - M[2]
        My = M[2] - M[0]
        Mz = M[0] - M[1]

        return Mx, My, Mz

    def check_convergence(self, est1, est2, tolerance):
        # checks if the algorithm has converged. est1 and est2 are two arrays
        difference =  abs(est2 - est1)
        print(np.sum(difference, axis=None))
        if np.sum(difference, axis=None) <= tolerance:
            return True
        else:
            return False

    def dist_to_charge(self, i, j):
        #calculates distance to center
        dx = i - self.dim / 2
        dy = j - self.dim / 2
        r = math.sqrt(dx**2 + dy**2)
        return r

    def Field_magnitude(self, field_x, field_y, field_z, i, j, k):
        # returns field magnitude
        return math.sqrt(field_x[i,j,k]**2 + field_y[i,j,k]**2 + field_z[i,j,k]**2)

    def get_results(self, kind):
        potential_data = []
        vector_data = []
        distance_to_charge = []

        if kind == 'E':
            field_x, field_y, field_z = self.Efield()
        elif kind == 'M':
            field_x, field_y, field_z = self.Mfield()

        for i in range(self.dim):
            for j in range(self.dim):
                for k in range(self.dim):
                    if k == int(self.dim / 2):
                        pot = [i, j, self.phi_latt[i,j,k]]
                        potential_data.append(pot)
                        magnitude = self.Field_magnitude(field_x, field_y, field_z, i, j, k)
                        vec = [i, j, field_x[i,j,k] / magnitude, field_y[i,j,k] / magnitude]
                        vector_data.append(vec)
                        dist = [self.dist_to_charge(i,j), self.phi_latt[i,j,k], magnitude]
                        distance_to_charge.append(dist)
        
        potential_data = np.array(potential_data)
        vector_data = np.array(vector_data)
        distance_to_charge = np.array(distance_to_charge)

        return potential_data, vector_data, distance_to_charge

    def plot_lattice(self, phi_sol, outdir, field):
        plt.figure(figsize=(8,8))
        plt.title(f"Solution to Poisson equation for {field} field")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.imshow(phi_sol[:][:][int(self.dim / 2)], cmap='gnuplot',)
        plt.colorbar()
        plt.savefig(f'{outdir}/phi_lattice_{self.dim}_{field}.png')
        plt.clf()

    def plot_quiver(self, vector_data, outdir, field):
        plt.figure(figsize=(8,8))
        plt.title(f"Vector Plot of Poisson equation for {field} field")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.quiver(vector_data[:,0], vector_data[:,1], vector_data[:,2], vector_data[:,3])
        plt.savefig(f'{outdir}/quiver_{self.dim}_{field}.png')

    def plot_distances(self, distance_to_charge, outdir, field_type):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (18,6))
        fig.suptitle(f'Behaviour of Potential and Field strength with distance to charge')

        x_log = np.log(distance_to_charge[:,0])
        pot_log = np.log(distance_to_charge[:,1])
        field_log = np.log(distance_to_charge[:,2])

        #concate all the logaritmic values of distance and potential
        cc = np.stack((x_log, pot_log))
        # finds indices for given range (100x100 lattice)
        indices = np.argwhere((cc[0] > 0.5) & (cc[0] < 2.5))
        #takes the corresponding values at indices found above
        reduced_cc = np.take(cc, indices, axis=1)
        # squeeze just reduces the dimension
        x_log_small = np.squeeze(reduced_cc[0,:,:])
        pot_log_small = np.squeeze(reduced_cc[1,:,:])

        # now similar for the field
        cc_f = np.stack((x_log, field_log))

        indices_f = np.argwhere((cc[0] > 2) & (cc[0] < 3)) 
        reduced_cc_f = np.take(cc_f, indices_f, axis=1)
        x_log_small_f = np.squeeze(reduced_cc_f[0,:,:])
        field_log_small = np.squeeze(reduced_cc_f[1,:,:])

        #create the linear fits
        p1 = np.polyfit(x_log_small, pot_log_small, 1)
        yfit1 = np.array(x_log_small)* p1[0] + p1[1]

        p2 = np.polyfit(x_log_small_f, field_log_small, 1)
        yfit2 = np.array(x_log_small_f)* p2[0] + p2[1]

        ax1.scatter(x_log, pot_log, s= 10)
        ax1.plot(x_log_small, yfit1, label = f"x(t) = {p1[0]:.3}t + {p1[1]:.3}", color='r')
        ax1.set_xlabel('Distance to charge (log scale)')
        ax1.set_ylabel('Potential (log scale)')
        ax1.legend()

        ax2.scatter(x_log, field_log, s= 10)
        ax2.plot(x_log_small_f, yfit2, label = f"x(t) = {p2[0]:.3}t + {p2[1]:.3}", color='r')
        ax2.set_xlabel('Distance to charge (log scale)')
        ax2.set_ylabel('Field (log scale)')
        ax2.legend()


        plt.tight_layout()
        fig.savefig(f"{outdir}/Potential_and_field_v_R_{self.dim}_{field_type}.png")


    @staticmethod
    def Solve_Poisson(lattice, algorithm, tolerance, kind, omega):

        convergence = False
        count = 0

        while (convergence == False):
            phi_old = lattice.phi_latt

            if algorithm == 'jacobi':
                lattice.phi_latt = lattice.jacobi()
            elif algorithm == 'gs':
                lattice.phi_latt = lattice.gauss_seidel(omega, SOR=False)
            elif algorithm == 'gs_SOR':
                lattice.phi_latt = lattice.gauss_seidel(omega, SOR=True)
                count += 1

            convergence = lattice.check_convergence(phi_old, lattice.phi_latt, tolerance)
        
        potential_data, vector_data, distance_to_charge = lattice.get_results(kind)

        return lattice.phi_latt, potential_data, vector_data, distance_to_charge, count

    
def main():

    args = parse_args()

    dim = args.dim 
    phi0 = args.phi0
    tolerance = args.tolerance
    omega = args.omega
    field_type = args.field_type
    SOR_calc = args.SOR_calc
    algorithm = args.algorithm


    if not os.path.isdir("Results"):
        os.makedirs("Results")
    outdir = 'Results'
    
    #search for optimum omega

    if SOR_calc:
        algorithm = "gs_SOR"
        omegas = np.arange(1.0, 2.0, 0.1)
        sweeps = np.empty(omegas.shape)

        for i, omega in enumerate(omegas):
            print('omega: {}'.format(omega))
            PS = poisson(dim, phi0)
            _ , _ , _ , _ , count = poisson.Solve_Poisson(PS, algorithm, tolerance, field_type, omega)
            sweeps[i] = count

        plt.figure(figsize=(10,6))
        plt.title(r'Optimum $\omega$ for SOR')
        plt.xlabel(r'$\omega$')
        plt.ylabel(r'Number of iterations')
        plt.plot(omegas, sweeps)
        plt.savefig(f'{outdir}/omegas_{dim}_{field_type}.png')

        with open(f"{outdir}/omegas_data_{dim}_{field_type}.dat", 'w+') as f:
            for i, j in zip(omegas, sweeps):
                f.write("{0:f} {1:f} \n".format(i, j))


    else:

        PS = poisson(dim, phi0)
        phi_sol, potential_data, vector_data, distance_to_charge, _ = poisson.Solve_Poisson(PS, algorithm, tolerance, field_type, omega)

        with open(f"{outdir}/Poisson_potential_data_{dim}_{algorithm}_{field_type}.dat", 'w+') as f:
            for i, j, phi_latt in potential_data:
                f.write("{0:f} {1:f} {2:f} \n".format(i, j, phi_latt))
        
        with open(f"{outdir}/Poisson_quiver_data_{dim}_{algorithm}_{field_type}.dat", 'w+') as f:
            for i, j, fieldx, fieldy in vector_data:
                f.write("{0:f} {1:f} {2:f} {3:f} \n".format(i, j, fieldx, fieldy))
                
        with open(f"{outdir}/Poisson_distance_data_{dim}_{algorithm}_{field_type}.dat", 'w+') as f:
            for dist, phi_latt, field in distance_to_charge:
                f.write("{0:f} {1:f} {2:f} \n".format(dist, phi_latt, field))


        PS.plot_lattice(phi_sol, outdir, field_type)
        PS.plot_quiver(vector_data, outdir, field_type)
        PS.plot_distances(distance_to_charge, outdir, field_type)


main()    

