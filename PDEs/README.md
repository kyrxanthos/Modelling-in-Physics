# Partial Differential Equations

## 1. **Cahn-Hilliard Equation**

### &diams; **Usage**

`python CH.py --help`

**Defaults**:   

- sweeps: 100000
- dim: 100
- kappa: 0.1
- a: 0.1
- M: 0.1
- delta_x: 1.5
- delta_y: 2.0
- phi0: 0.
- measurement animation



###  &diams; **Measurement**

The results are stored in a directory `Results\` where one can find the plots for the run and the .dat files which include all the measurements for each plot.



#### **Sample Results**

- Free Energy Plot phi0 = 0.0
> `FE_data_0.0.dat`, `FE_plot_dim_100_phi_0.0.png`

- Free Energy Plot phi0 = 0.5
> `FE_data_0.0.dat`, `FE_plot_dim_100_phi_0.0.png`

## 2. **Poisson Equation**

### &diams; **Usage**

- 3D: `python Poisson.py --help`

- 2D: `python Poisson2d.py --help`

**Defaults**:   

- dim: 50
- phi0: 0.
- tolerance: 0.01
- omega: 1.9
- field_type: 'E'
- SOR_calc: False
- algorithm: 'jacobi'

