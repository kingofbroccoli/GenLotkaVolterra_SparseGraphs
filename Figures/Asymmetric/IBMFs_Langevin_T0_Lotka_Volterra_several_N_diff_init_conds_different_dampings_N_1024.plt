set term postscript enhanced color eps dl 2.5
filenameoutput="IBMFs_Langevin_T0_Lotka_Volterra_several_N_diff_init_conds_different_dampings_N_1024.eps"
set output filenameoutput

reset
set size 1,1
set origin 0,0
set bmargin 3.2
set ylabel "{/=26 {/Symbol s}}" rotate by 90 offset 0,0
set xlabel "{/=26 {/Symbol m}}" offset 0,0.5
set key spacing 2.0 maxrows 5 width 8 at 0.29,0.13
set xra [-0.01:0.372]
set yra [-0.005:0.455]
set tics  font ",24"
set label "{/=26 SINGLE}" at 0.04, 0.27
set label "{/=26 FIXED}" at 0.047, 0.24
set label "{/=26 POINT}" at 0.046, 0.21


p "../../Data/Simulations/Asymmetric/Lotka-Volterra_transition_mult_epsilon_0.0_Partially_AsymGauss_lambda_1e-06_tol_1e-08_N_1024_c_3.00_T_0.0.txt" u ($1<0.354?$1:1/0):(($2+$3)/2):2:3 w yerror pt 6 lc 4 lw 2 title "{/=26 N=1024}" \
, "../../Data/IBMF/Asymmetric/IBMF_seq_RRG_T_0.000_lambda_0.000_PD_Lotka_Volterra_transitions_mult_av0_0.5_dn_0.5_ninitconds_10_tol_1e-6_maxiter_10000_eps_0.000_N_1024_c_3_damping_1.0_nseq_1.txt" u ($1<0.354?$1:1/0):(($2+$3)/2) w l lc 8 lw 2 dt 2 title "{/=26 IBMF d=1.0}" \
, "../../Data/IBMF/Asymmetric/IBMF_seq_RRG_T_0.000_lambda_0.000_PD_Lotka_Volterra_transitions_mult_av0_0.5_dn_0.5_ninitconds_10_tol_1e-6_maxiter_10000_eps_0.000_N_1024_c_3_damping_0.2_nseq_1.txt" u 1:(($2+$3)/2) w l lc 8 lw 3 dt 1  title "{/=26 IBMF d=0.2}"