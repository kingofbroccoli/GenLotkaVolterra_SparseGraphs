set term postscript enhanced color eps dl 2.5
filenameoutput="IBMF_seq_directed_ER_T0_Lotka_Volterra_several_mu_damping_02_with_inset.eps"
set output filenameoutput

reset
set multi
set size 1,1
set origin 0,0
set ylabel "{/=26 P_{nc}}" rotate by 90 offset 0,0
set xlabel "{/=26 c}" offset -3,0
set key spacing 2.0 maxrows 10 width 6 at 1.4,0.44
set xra [-0.002:4.0]
set yra [-0.01:1.01]
set tics  font ",24"
set label "{/=26 N = 65536}" at 2.2, 0.05


p "<awk '$1 ==65536 && $3==800 && $4==0' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_seq_directed_ER_summary_T_0.000_lambda_0.000_PD_Lotka_Volterra_final_av0_0.08_dn_0_ninitconds_1_tol_1e-6_maxiter_10000_eps_0.000_damping_0.2_nseq_1.txt" u ($2/1000):9:10 w l lc 1 lw 3 dt 1 title "{/=26 IBMF {/Symbol m} = 0.8}" \
, "<awk '$1 ==65536 && $3==900 && $4==0' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_seq_directed_ER_summary_T_0.000_lambda_0.000_PD_Lotka_Volterra_final_av0_0.08_dn_0_ninitconds_1_tol_1e-6_maxiter_10000_eps_0.000_damping_0.2_nseq_1.txt" u ($2/1000):9:10 w l lc 2 lw 3 dt 1 title "{/=26 {/Symbol m} = 0.9}" \
, "<awk '$1 ==65536 && $3==1100 && $4==0' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_seq_directed_ER_summary_T_0.000_lambda_0.000_PD_Lotka_Volterra_final_av0_0.08_dn_0_ninitconds_1_tol_1e-6_maxiter_10000_eps_0.000_damping_0.2_nseq_1.txt" u ($2/1000):9:10 w l lc 4 lw 3 dt 1 title "{/=26 {/Symbol m} = 1.1}" \
, "<awk '$1 ==65536 && $3==1500 && $4==0' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_seq_directed_ER_summary_T_0.000_lambda_0.000_PD_Lotka_Volterra_final_av0_0.08_dn_0_ninitconds_1_tol_1e-6_maxiter_10000_eps_0.000_damping_0.2_nseq_1.txt" u ($2/1000):9:10 w l lc 6 lw 3 dt 1 title "{/=26 {/Symbol m} = 1.5}" \
, "<awk '$1 ==65536 && $3==3000 && $4==0' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_seq_directed_ER_summary_T_0.000_lambda_0.000_PD_Lotka_Volterra_final_av0_0.08_dn_0_ninitconds_1_tol_1e-6_maxiter_10000_eps_0.000_damping_0.2_nseq_1.txt" u ($2/1000):9:10 w l lc 7 lw 3 dt 1 title "{/=26 {/Symbol m} = 3.0}" \
, 'p_fluc_directed_ER_mu_1.1.txt' w l lc 8 dt 2 lw 2 notitle \
, 'p_fluc_directed_ER_mu_1.5.txt' w l lc 8 dt 2 lw 2 notitle \
, 'p_fluc_directed_ER_mu_3.0.txt' w l lc 8 dt 2 lw 2 notitle


reset
set size 0.45,0.45
set origin 0.15,0.5
set ylabel "{/=22 P_{nc}}" rotate by 90 offset 0,0
set xlabel "{/=22 c}" offset -3,0
set key spacing 2.0 maxrows 10 width 6 at 1.3,0.85
set xra [-0.002:4.0]
set yra [-0.01:1.01]
set tics  font ",20"
set label "{/=22 N = 16384}" at 0.4, 0.8
set xtics 0,1,4

set arrow nohead from exp(1), -0.01 to exp(1),1.01 dt 2 lc 8 lw 2


p "<awk '$1 ==16384 && $3==800 && $4==0' ../../Data/Simulations/Directed_ErdosRenyi/Lotka-Volterra_summary_epsilon_0.000_Directed_Trivial_lambda_1e-06_tol_1e-08_T_0.000.txt" u ($2/1000):9:10 w yerror lc 1 lw 2 dt 3 pt 7 notitle \
, "<awk '$1 ==16384 && $3==900 && $4==0' ../../Data/Simulations/Directed_ErdosRenyi/Lotka-Volterra_summary_epsilon_0.000_Directed_Trivial_lambda_1e-06_tol_1e-08_T_0.000.txt" u ($2/1000):9:10 w yerror lc 2 lw 2 dt 3 pt 4 notitle \
, "<awk '$1 ==16384 && $3==1100 && $4==0' ../../Data/Simulations/Directed_ErdosRenyi/Lotka-Volterra_summary_epsilon_0.000_Directed_Trivial_lambda_1e-06_tol_1e-08_T_0.000.txt" u ($2/1000):9:10 w yerror lc 4 lw 2 dt 3 pt 8 notitle \
, "<awk '$1 ==16384 && $3==1500 && $4==0' ../../Data/Simulations/Directed_ErdosRenyi/Lotka-Volterra_summary_epsilon_0.000_Directed_Trivial_lambda_1e-06_tol_1e-08_T_0.000.txt" u ($2/1000):9:10 w yerror lc 6 lw 2 dt 3 pt 2 notitle \
, "<awk '$1 ==16384 && $3==3000 && $4==0' ../../Data/Simulations/Directed_ErdosRenyi/Lotka-Volterra_summary_epsilon_0.000_Directed_Trivial_lambda_1e-06_tol_1e-08_T_0.000.txt" u ($2/1000):9:10 w yerror lc 7 lw 2 dt 3 pt 10 notitle \
, "<awk '$1 ==16384 && $3==800 && $4==0' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_seq_directed_ER_summary_T_0.000_lambda_0.000_PD_Lotka_Volterra_final_av0_0.08_dn_0_ninitconds_1_tol_1e-6_maxiter_10000_eps_0.000_damping_0.2_nseq_1.txt" u ($2/1000):9 w l lc 1 lw 3 dt 1 notitle  \
, "<awk '$1 ==16384 && $3==900 && $4==0' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_seq_directed_ER_summary_T_0.000_lambda_0.000_PD_Lotka_Volterra_final_av0_0.08_dn_0_ninitconds_1_tol_1e-6_maxiter_10000_eps_0.000_damping_0.2_nseq_1.txt" u ($2/1000):9 w l lc 2 lw 3 dt 1 notitle \
, "<awk '$1 ==16384 && $3==1100 && $4==0' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_seq_directed_ER_summary_T_0.000_lambda_0.000_PD_Lotka_Volterra_final_av0_0.08_dn_0_ninitconds_1_tol_1e-6_maxiter_10000_eps_0.000_damping_0.2_nseq_1.txt" u ($2/1000):9 w l lc 4 lw 3 dt 1 notitle \
, "<awk '$1 ==16384 && $3==1500 && $4==0' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_seq_directed_ER_summary_T_0.000_lambda_0.000_PD_Lotka_Volterra_final_av0_0.08_dn_0_ninitconds_1_tol_1e-6_maxiter_10000_eps_0.000_damping_0.2_nseq_1.txt" u ($2/1000):9 w l lc 6 lw 3 dt 1 notitle \
, "<awk '$1 ==16384 && $3==3000 && $4==0' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_seq_directed_ER_summary_T_0.000_lambda_0.000_PD_Lotka_Volterra_final_av0_0.08_dn_0_ninitconds_1_tol_1e-6_maxiter_10000_eps_0.000_damping_0.2_nseq_1.txt" u ($2/1000):9 w l lc 7 lw 3 dt 1 notitle