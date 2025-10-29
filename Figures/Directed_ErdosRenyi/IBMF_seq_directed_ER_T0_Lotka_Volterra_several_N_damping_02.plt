set term postscript enhanced color eps dl 2.5
filenameoutput="IBMF_seq_directed_ER_T0_Lotka_Volterra_several_N_damping_02.eps"
set output filenameoutput

reset
set multi
set size 1,1
set origin 0,0
set ylabel "{/=26 P_{nc}}" rotate by 90 offset 0,0
set xlabel "{/=26 c}" offset -3,0
set key spacing 2.0 maxrows 10 width 6 at 2.1,0.8
set xra [1:4.0]
set yra [-0.01:1.01]
set tics  font ",24"
set arrow nohead from exp(1), -0.01 to exp(1),1.01 dt 2 lc 8


p "<awk '$1 ==8192 && $3==900 && $4==0' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_seq_directed_ER_summary_T_0.000_lambda_0.000_PD_Lotka_Volterra_final_av0_0.08_dn_0_ninitconds_1_tol_1e-6_maxiter_10000_eps_0.000_damping_0.2_nseq_1.txt" u ($2/1000):9:10 w yerrorl lc 7 lw 2 dt 1 pt 7 notitle \
, "<awk '$1 ==8192 && $3==1100 && $4==0' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_seq_directed_ER_summary_T_0.000_lambda_0.000_PD_Lotka_Volterra_final_av0_0.08_dn_0_ninitconds_1_tol_1e-6_maxiter_10000_eps_0.000_damping_0.2_nseq_1.txt" u ($2/1000):9:10 w yerrorl lc 7 lw 2 dt 3 pt 2 notitle \
, "dumb_data.txt" w l lc 7 lw 2 dt 1 title "{/=26 N=8192}" \
, "<awk '$1 ==16384 && $3==900 && $4==0' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_seq_directed_ER_summary_T_0.000_lambda_0.000_PD_Lotka_Volterra_final_av0_0.08_dn_0_ninitconds_1_tol_1e-6_maxiter_10000_eps_0.000_damping_0.2_nseq_1.txt" u ($2/1000):9:10 w yerrorl lc 2 lw 2 dt 1 pt 7 notitle \
, "<awk '$1 ==16384 && $3==1100 && $4==0' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_seq_directed_ER_summary_T_0.000_lambda_0.000_PD_Lotka_Volterra_final_av0_0.08_dn_0_ninitconds_1_tol_1e-6_maxiter_10000_eps_0.000_damping_0.2_nseq_1.txt" u ($2/1000):9:10 w yerrorl lc 2 lw 2 dt 3 pt 2 notitle \
, "dumb_data.txt" w l lc 2 lw 2 dt 1 title "{/=26 N=16384}" \
, "<awk '$1 ==32768 && $3==900 && $4==0' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_seq_directed_ER_summary_T_0.000_lambda_0.000_PD_Lotka_Volterra_final_av0_0.08_dn_0_ninitconds_1_tol_1e-6_maxiter_10000_eps_0.000_damping_0.2_nseq_1.txt" u ($2/1000):9:10 w yerrorl lc 4 lw 2 dt 1 pt 7 notitle \
, "<awk '$1 ==32768 && $3==1100 && $4==0' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_seq_directed_ER_summary_T_0.000_lambda_0.000_PD_Lotka_Volterra_final_av0_0.08_dn_0_ninitconds_1_tol_1e-6_maxiter_10000_eps_0.000_damping_0.2_nseq_1.txt" u ($2/1000):9:10 w yerrorl lc 4 lw 2 dt 3 pt 2 notitle \
, "dumb_data.txt" w l lc 4 lw 2 dt 1 title "{/=26 N=32768}" \
, "<awk '$1 ==65536 && $3==900 && $4==0' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_seq_directed_ER_summary_T_0.000_lambda_0.000_PD_Lotka_Volterra_final_av0_0.08_dn_0_ninitconds_1_tol_1e-6_maxiter_10000_eps_0.000_damping_0.2_nseq_1.txt" u ($2/1000):9:10 w yerrorl lc 6 lw 2 dt 1 pt 7 notitle \
, "<awk '$1 ==65536 && $3==1100 && $4==0' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_seq_directed_ER_summary_T_0.000_lambda_0.000_PD_Lotka_Volterra_final_av0_0.08_dn_0_ninitconds_1_tol_1e-6_maxiter_10000_eps_0.000_damping_0.2_nseq_1.txt" u ($2/1000):9:10 w yerrorl lc 6 lw 2 dt 3 pt 2 notitle \
, "dumb_data.txt" w l lc 6 lw 2 dt 1 title "{/=26 N=65536}" \
, 'dumb_data.txt' w yerrorl lc 8 lw 2 dt 1 pt 7 title "{/=26 {/Symbol m} = 0.9}" \
, 'dumb_data.txt' w yerrorl lc 8 lw 2 dt 3 pt 2 title "{/=26 {/Symbol m} = 1.1}"# \
#, 'p_fluc_directed_ER_mu_1.1.txt' w l lc 8 dt 2 notitle