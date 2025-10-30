# Figures/Symmetric — gnuplot scripts and data dependencies

This folder contains gnuplot scripts (.plt) that generate EPS figures for the symmetric-panel figures. The scripts read data from `Data/BP/Symmetric` and `Data/IBMF/Symmetric` (and a few local helper files). Run gnuplot from the repository root so the relative paths in the scripts resolve.

Prerequisites
- gnuplot (5.x recommended)
- POSIX shell (Linux)

Run
- Single script (from repo root):
```bash
gnuplot Figures/Symmetric/<script>.plt
```
- All scripts:
```bash
for f in Figures/Symmetric/*.plt; do
  echo "Running $f"
  gnuplot "$f"
done
```

Inspect data files (headers / sample rows)
```bash
head -n 20 Data/BP/Symmetric/*.txt
head -n 20 Data/IBMF/Symmetric/*.txt
```
If header lines are commented, use:
```bash
awk '/^[^#]/ {print; exit}' <filename>
```

Per-script summary (files referenced — inferred)

1) BP_sigma0_Gaussianity_T_003_mu_004.plt
- Output: BP_sigma0_Gaussianity_T_003_mu_004.eps
- Data files referenced:
  - ../../Data/BP/Symmetric/BP_cont_LV_sigma0_true_marginals_gaussian_part_N_128_T_0.030_mu_0.040.txt
  - ../../Data/BP/Symmetric/BP_cont_LV_sigma0_true_marginals_fit_trunc_gauss_N_128_T_0.030_mu_0.040.txt
  - ../../Data/BP/Symmetric/BP_cont_LV_sigma0_diff_true_marginals_fit_trunc_gauss_N_128_T_0.030_mu_0.040.txt

2) BP_sigma0_Gaussianity_T_003_mu_006.plt
- Output: BP_sigma0_Gaussianity_T_003_mu_006.eps
- Data files referenced:
  - ../../Data/BP/Symmetric/BP_cont_LV_sigma0_true_marginals_gaussian_part_N_128_T_0.030_mu_0.060.txt
  - ../../Data/BP/Symmetric/BP_cont_LV_sigma0_true_marginals_fit_trunc_gauss_N_128_T_0.030_mu_0.060.txt
  - ../../Data/BP/Symmetric/BP_cont_LV_sigma0_diff_true_marginals_fit_trunc_gauss_N_128_T_0.030_mu_0.060.txt

3) BP_sigma0_Gaussianity_T_003_mu_012.plt
- Output: BP_sigma0_Gaussianity_T_003_mu_012.eps
- Data files referenced:
  - ../../Data/BP/Symmetric/BP_cont_LV_sigma0_true_marginals_gaussian_part_N_128_T_0.030_mu_0.120.txt
  - ../../Data/BP/Symmetric/BP_cont_LV_sigma0_true_marginals_fit_trunc_gauss_N_128_T_0.030_mu_0.120.txt
  - ../../Data/BP/Symmetric/BP_cont_LV_sigma0_diff_true_marginals_fit_trunc_gauss_N_128_T_0.030_mu_0.120.txt

4) IBMFs_BP_sigma0_from_input_N_1024.plt
- Output: IBMFs_BP_sigma0_from_input_N_1024.eps
- Data files referenced:
  - ../../Data/IBMF/Symmetric/IBMF_seq_from_input_RRG_PD_Lotka_Volterra_summary_av0_0.08_lambda_1e-6_tol_1e-6_maxiter_10000_eps_1.0_sigma_0.000_N_1024_c_3_damping_1.0_nseq_10000.txt
  - ../../Data/BP/Symmetric/BP_from_input_RRG_PD_Lotka_Volterra_summary_N_1024_damping_1.0.txt

5) IBMFs_BP_sigma0_Lotka_Volterra_changing_T.plt
- Output: IBMFs_BP_sigma0_Lotka_Volterra_changing_T.eps
- Data files referenced:
  - ../../Data/IBMF/Symmetric/IBMF_seq_RRG_PD_Lotka_Volterra_transitions_av0_0.08_lambda_1e-6_tol_1e-6_maxiter_10000_eps_1.000_sigma_0.000_N_1024_c_3_damping_1.0_nseq_10.txt
  - ../../Data/BP/Symmetric/BP_transition_sigma0_N_128_sequential_changing_order.txt
  - Stars_In_IBMFs_BP_sigma0.txt (local helper file)

Helper files in this folder
- Stars_In_IBMFs_BP_sigma0.txt — local helper with plotting markers/positions used by IBMFs_BP_sigma0_Lotka_Volterra_changing_T.plt
