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
- `seed`: Optional random seed

**Output:**
- `BP_cont_LV_sigma0_true_marginals_N_*_T_*_mu_*.txt`: Node marginal distributions
- `BP_cont_LV_sigma0_true_marginals_gaussian_part_N_*_T_*_mu_*.txt`: Gaussian component analysis

#### [`cont-BP-grafo-transition-line.c`](cont-BP-grafo-transition-line.c)
Implements continuous BP to detect phase transitions by scanning parameter space (temperature T and interaction strength MU).

**Key Features:**
- Scans temperature and interaction strength space
- Detects convergence for each parameter combination
- Identifies critical points where BP fails to converge

**Usage:**
```bash
./program MU T N_tot [seed]
```

**Output:**
- `MU_crit.dat`: Critical points where BP fails to converge

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

**Mathematical Details:**
- Uses parabolic cylinder functions for initial BP message
- Implements truncated Gaussian approximation
- Grid spacing: $\Delta = \frac{\text{step}}{50}$ for fine resolution

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

## Algorithm Overview

### Belief Propagation on Random Regular Graphs

The code implements continuous BP for a factor graph defined on random 3-regular graphs. Key equations:

**Message Update:**
$$\eta_{i \to j}(n_i) \propto e^{-\beta(n_i^2 - 2n_i)} \prod_{k \in \partial i \setminus j} \eta_{k \to i}(n_i)$$

**Marginal Distribution:**
$$\eta_i(n_i) \propto e^{-\beta(n_i^2 - 2n_i)} \prod_{k \in \partial i} \eta_{k \to i}(n_i)$$

**Numerical Integration:**
Uses Simpson's rule with adaptive grid discretization.

## VS Code Configuration

The repository includes VS Code settings:
- `.vscode/launch.json`: LLDB debugging configuration
- `.vscode/c_cpp_properties.json`: Include paths and compiler settings
- `.vscode/settings.json`: C++ runner and warning configurations

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

## Output Formats

**Distribution Files:**
```
# n_i    av(eta)    eta_i[0]    eta_i[1]    ...    eta_i[N-1]
0.0      0.001234   0.001200    0.001150    ...    0.001300
0.002    0.001256   0.001230    0.001180    ...    0.001320
...
```

**Critical Points File (`MU_crit.dat`):**
```
T_value    MU_critical
0.001      0.123
0.002      0.145
...
```

## References

The implementation is based on Belief Propagation algorithms for continuous variables on factor graphs and phase transition analysis in random graph models.
