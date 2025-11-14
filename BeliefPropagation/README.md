# Belief Propagation on Random Regular Graphs - Continuous Analysis

This repository contains C and Python implementations for analyzing continuous Belief Propagation (BP) on random regular graphs, with focus on distribution analysis and phase transition detection.

## Files Overview

### C Programs

#### [`cont-BP-grafo_print_distr.c`](cont-BP-grafo_print_distr.c)
Continuous Belief Propagation implementation that computes and outputs marginal distributions for nodes in a random regular graph.

**Key Features:**
- Implements continuous BP algorithm with message passing
- Computes true marginal distributions for each node
- Separates analysis into total distribution and Gaussian part
- Uses Simpson's rule for numerical integration
- Outputs detailed statistics: mean, variance, skewness, and kurtosis

**Usage:**
```bash
./program MU T N_tot N damping [seed]
```
- `MU`: Interaction strength parameter
- `T`: Temperature parameter
- `N_tot`: Grid size for discretization
- `N`: Number of nodes
- `damping`: Damping factor for message updates (0-1)
- `seed` (optional): Random seed for graph generation

**Output:**
- `BP_cont_LV_sigma0_true_marginals_N_*_T_*_mu_*.txt`: Node marginal distributions
- `BP_cont_LV_sigma0_true_marginals_gaussian_part_N_*_T_*_mu_*.txt`: Gaussian component analysis

#### [`cont-BP-grafo-transition-line.c`](cont-BP-grafo-transition-line.c)
Implements continuous Belief Propagation on a random regular graph to detect multiple-equilibria phase transitions by scanning (T, MU) parameter space.

**Key Features:**
- Scans temperature and interaction strength space
- Detects convergence for each parameter combination
- Identifies critical points where BP fails to converge

**Usage:**
```bash
./program MU_start T_start N_tot [seed]
```
- `MU_start`: Starting value of the interaction parameter MU. The program increases MU in steps of 0.001.
- `T_start`: Starting temperature. The program increases T until 0.0342, using a step that is 0.001 or 0.0005 depending on T.
- `N_tot`: Grid size for discretization
- `seed` (optional): Random seed for graph generation

**Output:**
- `MU_crit.dat`: Critical points (T, MU) where BP fails to converge

#### [`functions-cont-BP-grafo-transition-line.h`](functions-cont-BP-grafo-transition-line.h)
Header file containing core BP functions shared between C programs.

**Main Functions:**
- `allocate_tree()`, `allocate_conf()`, `allocate_par()`, `allocate_obs()`: Memory allocation
- `create_grid()`: Discretizes the domain $[0, 1.2]$
- `compute_simpson_weights()`: Computes weights for Simpson's rule integration
- `compute_field()`: Computes external field
- `compute_exp_interaction()`: Computes interaction terms
- `BP_iteration()`: Single BP update step
- `compute_eta_hat_MOMENTS()`: Computes moments from marginals
- `MakeRandomRegularGraph()`: Generates random 3-regular graphs

### Python Script

#### [`corrected_skewness.py`](corrected_skewness.py)
Post-processing tool to fit truncated Gaussian distributions to BP results and analyze higher moments.

**Key Features:**
- Reads BP output distributions
- Fits truncated Gaussian models
- Computes skewness and kurtosis corrections
- Generates comparison plots and error analysis
- Uses `scipy.stats.truncnorm` for truncated Gaussian fitting

**Main Functions:**
- `read_distribution()`: Parse output files
- `compute_moments()`: Calculate mean, variance, skewness, kurtosis via Simpson integration
- `fit_truncated_gaussian()`: Fit parameters using `fsolve`
- `plot_distribution()`: Visualize BP vs. fitted Gaussian
- `plot_difference()`: Show relative errors

**Usage:**
```python
python corrected_skewness.py
```
(Configure paths and parameters in `main()`)

## Requirements

### C Programs
- GCC/Clang compiler
- GSL (GNU Scientific Library) for special functions:
  - `gsl_sf_gamma()`: Gamma function
  - `gsl_sf_hyperg_1F1()`: Hypergeometric function

**Compile:**
```bash
gcc -o program program.c -lgsl -lgslcblas -lm
```

### Python Script
- Python 3.7+
- NumPy
- SciPy
- Matplotlib

**Install dependencies:**
```bash
pip install numpy scipy matplotlib
```

## Parameters

| Parameter | Range | Description |
|-----------|-------|-------------|
| `MU` (`μ`) | 0.001-0.5 | Coupling strength |
| `T` | 0.001-0.05 | Temperature (inverse: $\beta = 1/T$) |
| `N` | Typical: 128, 256 | Number of graph nodes |
| `N_tot` | 1000-2000 | Grid discretization points |
| `damping` | 0-1 | Message damping factor |
| `LAMBDA` | Fixed at $10^{-6}$ | Power law parameter |
| `C` | Fixed at 3 | Graph connectivity (regular 3-graphs) |


Authors / Contact
-----------------
This code is part of the `GenLotkaVolterra_SparseGraphs` repository.
If you need help or wish to improve documentation, open an issue or contact the maintainer.

License
-------
Check the repository root for license information; this folder does not contain a separate license file.
