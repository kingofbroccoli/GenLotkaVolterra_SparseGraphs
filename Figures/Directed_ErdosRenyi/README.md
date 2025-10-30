# Figures/Directed_ErdosRenyi — gnuplot scripts and data dependencies

This folder contains gnuplot scripts (.plt) that generate EPS figures for the Section 4.B of the paper (directed Erdős–Rényi). The scripts read data from `Data/IBMF/Directed_ErdosRenyi` and `Data/Simulations/Directed_ErdosRenyi`. Run gnuplot from the repository root so the relative paths in the scripts resolve.

Prerequisites
- gnuplot (5.x recommended)
- POSIX shell (Linux)

Run
- Single script (from repo root):
```bash
gnuplot Figures/Directed_ErdosRenyi/<script>.plt
```
- All scripts:
```bash
for f in Figures/Directed_ErdosRenyi/*.plt; do
  echo "Running $f"
  gnuplot "$f"
done
```

Inspect data files (headers / sample rows)
```bash
head -n 20 Data/IBMF/Directed_ErdosRenyi/*.txt
head -n 20 Data/Simulations/Directed_ErdosRenyi/*.txt
```
If header lines are commented, use:
```bash
awk '/^[^#]/ {print; exit}' <filename>
```

Per-script summary (files referenced — inferred)

1) IBMFs_vs_Langevin_directed_ER_T0_Lotka_Volterra_N_8192_runtimes_damping_02.plt
- Output: IBMFs_vs_Langevin_directed_ER_T0_Lotka_Volterra_N_8192_runtimes_damping_02.eps
- Data files referenced:
  - ../../Data/Simulations/Directed_ErdosRenyi/Lotka-Volterra_summary_epsilon_0.000_Directed_Trivial_lambda_1e-06_tol_1e-08_T_0.000.txt
  - ../../Data/IBMF/Directed_ErdosRenyi/IBMF_seq_directed_ER_summary_T_0.000_lambda_0.000_PD_Lotka_Volterra_final_av0_0.08_dn_0_ninitconds_1_tol_1e-6_maxiter_10000_eps_0.000_damping_0.2_nseq_1_Dariah.txt

2) IBMF_seq_directed_ER_T0_Lotka_Volterra_several_sigma_damping_02.plt
- Output: IBMF_seq_directed_ER_T0_Lotka_Volterra_several_sigma_damping_02.eps
- Data files referenced:
  - ../../Data/IBMF/Directed_ErdosRenyi/IBMF_seq_directed_ER_summary_*.txt

3) IBMF_seq_directed_ER_T0_Lotka_Volterra_several_N_damping_02.plt
- Output: IBMF_seq_directed_ER_T0_Lotka_Volterra_several_N_damping_02.eps
- Data files referenced:
  - ../../Data/IBMF/Directed_ErdosRenyi/IBMF_seq_directed_ER_summary_*.txt
  - ../../Data/Simulations/Directed_ErdosRenyi/Lotka-Volterra_summary_epsilon_0.000_Directed_Trivial_*.txt

4) IBMF_seq_directed_ER_T0_Lotka_Volterra_several_mu_damping_02_with_inset.plt
- Output: IBMF_seq_directed_ER_T0_Lotka_Volterra_several_mu_damping_02_with_inset.eps
- Data files referenced:
  - ../../Data/IBMF/Directed_ErdosRenyi/IBMF_seq_directed_ER_summary_*.txt
  - ../../Data/Simulations/Directed_ErdosRenyi/Lotka-Volterra_summary_epsilon_0.000_Directed_Trivial_*.txt (used for inset)
  - Figures/Directed_ErdosRenyi/p_fluc_directed_ER_mu_1.1.txt
  - Figures/Directed_ErdosRenyi/p_fluc_directed_ER_mu_1.5.txt
  - Figures/Directed_ErdosRenyi/p_fluc_directed_ER_mu_3.0.txt

5) IBMF_directed_ER_T0_Lotka_Volterra_several_N.plt
- Output: IBMF_directed_ER_T0_Lotka_Volterra_several_N.eps
- Data files referenced:
  - ../../Data/IBMF/Directed_ErdosRenyi/IBMF_T0_directed_ER_PD_Lotka_Volterra_summary_av0_0.08_tol_1e-4_maxiter_10000.txt
  - Figures/Directed_ErdosRenyi/p_fluc_IBMF_directed_ER_n0_2.txt

6) IBMF_directed_ER_T0_Lotka_Volterra_N_65536_several_mu.plt
- Output: IBMF_directed_ER_T0_Lotka_Volterra_N_65536_several_mu.eps
- Data files referenced:
  - ../../Data/IBMF/Directed_ErdosRenyi/IBMF_T0_directed_ER_PD_Lotka_Volterra_summary_av0_0.08_tol_1e-4_maxiter_10000.txt
  - Figures/Directed_ErdosRenyi/p_fluc_IBMF_directed_ER_n0_2.txt