# Figures/Asymmetric - Gnuplot Scripts Documentation

This folder contains gnuplot scripts (.plt) that generate EPS figures for the Section 4.A (with asymmetric graphs) of the paper. Scripts use relative paths to data files in the `Data/` directory.

## Prerequisites
- gnuplot (5.x recommended)
- POSIX shell (Linux)

## Usage
Run from repository root:
```bash
# Run single script
gnuplot Figures/Asymmetric/<script>.plt

# Run all scripts in folder
for f in Figures/Asymmetric/*.plt; do
  echo "Running $f"
  gnuplot "$f"
done
```

## Scripts and Their Data Dependencies

### Langevin_T0_Lotka_Volterra_crossover.plt
- Output: Langevin_T0_Lotka_Volterra_crossover.eps
- Data files:
  ```
  Data/Simulations/Asymmetric/Lotka-Volterra_transition_mult_epsilon_0.0_Partially_AsymGauss_lambda_1e-06_tol_1e-08_N_1024_c_3.00_T_0.0.txt
  Data/Simulations/Asymmetric/Lotka-Volterra_transition_div_epsilon_0.0_Partially_AsymGauss_lambda_1e-06_tol_1e-08_N_1024_c_3.00_T_0.0.txt
  ```
- Column usage:
  - col1: x (parameter μ)
  - col2: lower transition bound
  - col3: upper transition bound
  - Plotted as: y = (col2+col3)/2 with error bars from col2 to col3

### IBMFs_Langevin_T0_Lotka_Volterra_several_N_divergence.plt
- Output: IBMFs_Langevin_T0_Lotka_Volterra_several_N_divergence.eps
- Data files:
  - Figures/Asymmetric/dumb_data.txt (reference curve)
  - Transition files for N = 256,512,1024,2048,4096:
    ```
    Data/Simulations/Asymmetric/Lotka-Volterra_transition_div_*_N_<N>.txt
    Data/IBMF/Asymmetric/IBMF_T0_seq_RRG_PD_Lotka_Volterra_transitions_div_*_N_<N>.txt
    ```

### IBMFs_Langevin_T0_Lotka_Volterra_several_N_divergence_allrange.plt
- Output: IBMFs_Langevin_T0_Lotka_Volterra_several_N_divergence_allrange.eps 
- Uses same data files as above script
- Additional feature: Vertical reference line at x = -1/3

### IBMFs_Langevin_T0_Lotka_Volterra_several_N_diff_init_conds.plt
- Output: IBMFs_Langevin_T0_Lotka_Volterra_several_N_diff_init_conds.eps
- Data files:
  - Figures/Asymmetric/dumb_data.txt
  - For N = 256,512,1024,2048,4096:
    ```
    Data/Simulations/Asymmetric/Lotka-Volterra_transition_mult_*_N_<N>.txt
    Data/IBMF/Asymmetric/IBMF_seq_RRG_T_*_transitions_mult_*_N_<N>.txt
    ```
- Column usage:
  - Transition files: col1=x, col2=lower, col3=upper
  - Some files include col4 (multiplicity flag)

### IBMFs_Langevin_Lotka_Volterra_abundances_single_graph.plt
- Output: IBMFs_Langevin_Lotka_Volterra_abundances_single_graph.eps
- Data files:
  - Figures/Asymmetric/IBMF_vs_Langevin_single_abundances.txt
    - col1: n^{SIM}
    - col4: n^{IBMF}
  - Time series:
    ```
    Data/Simulations/Asymmetric/One_graph/Lotka-Volterra_Extraction_1_Measure_2_*.txt
    ```
    - col1: time
    - col2+: species abundances

## Helper Files
- IBMF_vs_Langevin_single_abundances.txt: Abundance comparison data