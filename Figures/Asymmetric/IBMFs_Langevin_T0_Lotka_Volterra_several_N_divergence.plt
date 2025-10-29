set term postscript enhanced color eps dl 2.5
filenameoutput="IBMFs_Langevin_T0_Lotka_Volterra_several_N_divergence.eps"
set output filenameoutput

reset
set size 1,1
set origin 0,0
set ylabel "{/=26 {/Symbol s}}" rotate by 90 offset 0,0
set xlabel "{/=26 {/Symbol m}}" offset 0,0.5
set key spacing 2.0 maxrows 5 width 5 at 0.37,0.44
set xra [-0.01:0.372]
set yra [0.35:0.71]
set tics  font ",24"
set label "{/=26 UNBOUNDED}" at 0.05, 0.65
set label "{/=26 GROWTH}" at 0.065, 0.625


p "../../Data/Simulations/Asymmetric/Lotka-Volterra_transition_div_epsilon_0.0_Partially_AsymGauss_lambda_1e-06_tol_1e-08_N_4096_c_3.00_T_0.0.txt" u ($1>0?$1:1/0):(($2+$3)/2):2:3 w yerror pt 2 lc 3 lw 2 title "{/=26 N=4096}" \
, "../../Data/Simulations/Asymmetric/Lotka-Volterra_transition_div_epsilon_0.0_Partially_AsymGauss_lambda_1e-06_tol_1e-08_N_2048_c_3.00_T_0.0.txt" u ($1>0?$1:1/0):(($2+$3)/2):2:3 w yerror pt 4 lc 7 lw 2 title "{/=26 N=2048}" \
, "../../Data/Simulations/Asymmetric/Lotka-Volterra_transition_div_epsilon_0.0_Partially_AsymGauss_lambda_1e-06_tol_1e-08_N_1024_c_3.00_T_0.0.txt" u ($1>0?$1:1/0):(($2+$3)/2):2:3 w yerror pt 6 lc 4 lw 2 title "{/=26 N=1024}" \
, "../../Data/Simulations/Asymmetric/Lotka-Volterra_transition_div_epsilon_0.0_Partially_AsymGauss_lambda_1e-06_tol_1e-08_N_512_c_3.00_T_0.0.txt" u ($1>0?$1:1/0):(($2+$3)/2):2:3 w yerror pt 8 lc 2 lw 2 title "{/=26 N=512}" \
, "../../Data/Simulations/Asymmetric/Lotka-Volterra_transition_div_epsilon_0.0_Partially_AsymGauss_lambda_1e-06_tol_1e-08_N_256_c_3.00_T_0.0.txt" u ($1>0?$1:1/0):(($2+$3)/2):2:3 w yerror pt 10 lc 1 lw 2 title "{/=26 N=256}" \
, "../../Data/IBMF/Asymmetric/IBMF_T0_seq_RRG_PD_Lotka_Volterra_transitions_div_av0_0.08_tol_1e-6_maxiter_10000_eps_0.000_N_4096_c_3_damping_0.2_nseq_10.txt" u ($1>0?$1:1/0):(($2+$3)/2) w l lc 8 lw 3 dt 1 notitle \
, "../../Data/IBMF/Asymmetric/IBMF_T0_seq_RRG_PD_Lotka_Volterra_transitions_div_av0_0.08_tol_1e-6_maxiter_10000_eps_0.000_N_2048_c_3_damping_0.2_nseq_10.txt" u ($1>0?$1:1/0):(($2+$3)/2) w l lc 8 lw 3 dt 1 notitle \
, "../../Data/IBMF/Asymmetric/IBMF_T0_seq_RRG_PD_Lotka_Volterra_transitions_div_av0_0.08_tol_1e-6_maxiter_10000_eps_0.000_N_1024_c_3_damping_0.2_nseq_10.txt" u ($1>0?$1:1/0):(($2+$3)/2) w l lc 8 lw 3 dt 1 notitle \
, "../../Data/IBMF/Asymmetric/IBMF_T0_seq_RRG_PD_Lotka_Volterra_transitions_div_av0_0.08_tol_1e-6_maxiter_10000_eps_0.000_N_512_c_3_damping_0.2_nseq_10.txt" u ($1>0?$1:1/0):(($2+$3)/2) w l lc 8 lw 3 dt 1 notitle \
, "../../Data/IBMF/Asymmetric/IBMF_T0_seq_RRG_PD_Lotka_Volterra_transitions_div_av0_0.08_tol_1e-6_maxiter_10000_eps_0.000_N_256_c_3_damping_0.2_nseq_10.txt" u ($1>0?$1:1/0):(($2+$3)/2) w l lc 8 lw 3 dt 1 notitle \
, "dumb_data.txt" u 1:2 w l lc 8 lw 3 dt 1 title "{/=26 IBMF}"