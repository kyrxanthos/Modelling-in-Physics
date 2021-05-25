# Ising Model

## &diams; Usage

`python ising.py --help`

**Defaults**:   

    -dynamics Glauber
    -sweeps 10000
    -dim 50
    -animation False
    -temp 1.2

### **To show animation**

`python ising.py --animate True`

**Example**:

`python ising.py --animate True --dynamics Kawasaki --sweeps 1000 --temp 1.5`

### **For measurements**

**Example**:

`python ising.py  --dynamics Glauber --sweeps 1000 --dim 30`

**Note**: The program is coded in such a way that does not allow the user to both show animation and make measurements. This was intentional since the animation takes a lot of time to finish and measurements are faster without the animation.

## &diams; Measurement Results

The results are stored in a directory `Results\` where you can find the plots for the run and the .dat files which include all the measurements for each temperature point. 

### **.dat files sequence**

Temperature, Energy Mean, Specific Heat Capacity , Magnetisation Mean, Susceptibility, Error of Susceptibility, Error of Heat Capacity

### **Sample Results**

- Glauber Dynamics, 10,000 sweeps, lattice 50x50 
> `measurements_Glauber_10000_50.dat`, `plots_Glauber_10000.png`

- Kawasaki Dynamics, 10,000 sweeps, lattice 50x50
> `measurements_Kawasaki_10000_50`, `plots_Kawasaki_10000.png`