# Celular Automata

## 1. **Game of Life**

### &diams; **Usage**

`python GOL.py --help`

**Defaults**:   

    -init random
    -sweeps 10000
    -dim 50
    -measurement animation

#### **Animation**

- `python GOL.py --measurement animation --init random`
- `python GOL.py --measurement animation --init glider`
- `python GOL.py --measurement animation --init oscilator`
- `python GOL.py --measurement animation --init absorbing`

#### **Equilibration**

`python GOL.py --measurement equilibration --sweeps  100`


#### **Glider Speed**

`python GOL.py --measurement glider_speed --sweeps  1000`


###  &diams; **Measurement**

The results are stored in a directory `Results\` where one can find the plots for the run and the .dat files which include all the measurements for each plot.

**Datafiles Key:**

- equilibration_data.dat columns:
> `time`

- histogram_data_dim_50.dat columns:
> `Bin, Frequency`

- equilibration_data_dim_50.dat columns:
> `times, x position, y position, r position`

**Note**: The program is coded in such a way that does not allow the user to both show animation and make measurements. This was intentional since the animation takes a lot of time to finish and measurements are faster without the animation.


#### **Example Results**

- Glider speed, 1,000 sweeps, lattice 50x50 
> `glider_speed_data.dat`, `glider_speed_plots_dim_50.png`

- Equilibration histogram, 300 sweeps, lattice 50x50
> `histogram_data_dim_50.dat`, `equilibration_data_dim_50.dat` `histogram_dim_50.png`

## 2. **SIRS**

### &diams; **Usage**

`python SIRS.py --help`

**Defaults**:   

    -sweeps 10000
    -dim 50
    -measurement animation

#### **Animation**

- Absorbing state <br />
`python SIRS.py --dim 50 --measurement animation --p1  0.1 --p2 0.5 --p3 0.01`

- Dynamic Equilibrium <br />
`python SIRS.py --dim 50 --measurement animation --p1  0.5 --p2 0.5 --p3 0.5`

- Cyclic Wave of Infections <br />
`python SIRS.py --dim 100 --measurement animation --p1  0.8 --p2 0.1 --p3 0.01`

#### **Phase Diagram**

`python SIRS.py --dim 50 --measurement phase_space --sweeps 1000`


#### **Variance with fixed p3**

`python SIRS.py --dim 50 --measurement variance --sweeps 10000`

#### **Vaccine**

- Vaccination with p1 = p2 = p3 = 0.5: <br />
`python SIRS.py --dim 50 --measurement vaccine --p1 0.5 --p2 0.5 --p3 0.5`

- Vaccination with p1 = 0.8 p2 = 0.1 p3 = 0.02: <br />
`python SIRS.py --dim 100 --measurement vaccine --p1 0.8 --p2 0.1 --p3 0.02`

###  &diams; **Measurement**

The results are stored in a directory `Results\` where one can find the plots for the run and the .dat files which include all the measurements for each plot.

**Datafiles Key:**

- phase_space_data.dat columns:
> `p1, p3, Mean(<I>), Var(<I>)`

- variance_data.dat columns:
> `p1, Var(<I>), error`

- vaccine_data.dat columns:
> `f_Im, Mean(<I>), error`

**Note**: The program is coded in such a way that does not allow the user to both show animation and make measurements. This was intentional since the animation takes a lot of time to finish and measurements are faster without the animation.


#### **Sample Results**

- Phase Space (Mean), 1,000 sweeps, lattice 50x50 
> `phase_space_data.dat`, `phase_diagram_Mean.png`

- Phase Space (Variance), 1,000 sweeps, lattice 50x50 
> `phase_space_data.dat`, `phase_diagram_Variance.png`

- Variance Cut, 10,000 sweeps, lattice 50x50 
> `variance_data.dat`, `variance_plot_dim50_p3_0.5.png`

- Immunity plots, 1,000 sweeps, lattice 50x50, p1 = p2 = p3 = 0.5
> `vaccine_dim50_data.dat`, `vaccine_plot_dim50_p3_0.5.png`

- Immunity plots, 1,000 sweeps, lattice 50x50, p1 = 0.8 p2 = 0.1 p3 = 0.02
> `vaccine_dim100_data.dat`, `vaccine_plot_dim100_p3_0.02.png`
