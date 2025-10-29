set term postscript enhanced color eps dl 2.5
filenameoutput="IBMF_seq_directed_ER_T0_Lotka_Volterra_several_sigma_damping_02.eps"
set output filenameoutput

reset
set multi
set size 1,1
set origin 0,0
set ylabel "{/=26 P_{nc}}" rotate by 90 offset 0,0
set xlabel "{/=26 c}" offset 0,0
set key spacing 2.0 maxrows 10 width 6 at 4.88,0.8
set xra [4.5:5.7]
set yra [-0.01:1.01]
set tics  font ",24"

p "<awk '$1 ==8192 && $3==700 && $4==0' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_seq_directed_ER_summary_T_0.000_lambda_0.000_PD_Lotka_Volterra_final_av0_0.08_dn_0_ninitconds_1_tol_1e-6_maxiter_10000_eps_0.000_damping_0.2_nseq_1.txt" u ($2/1000):9:10 w yerrorl lc 7 lw 2 dt 1 pt 7 notitle \
, "<awk '$1 ==8192 && $3==700 && $4==150' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_seq_directed_ER_summary_T_0.000_lambda_0.000_PD_Lotka_Volterra_final_av0_0.08_dn_0_ninitconds_1_tol_1e-6_maxiter_10000_eps_0.000_damping_0.2_nseq_1.txt" u ($2/1000):9:10 w yerrorl lc 7 lw 2 dt 3 pt 2 notitle \
, "dumb_data.txt" w l lc 7 lw 2 dt 1 title "{/=26 N=8192}" \
, "<awk '$1 ==16384 && $3==700 && $4==0' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_seq_directed_ER_summary_T_0.000_lambda_0.000_PD_Lotka_Volterra_final_av0_0.08_dn_0_ninitconds_1_tol_1e-6_maxiter_10000_eps_0.000_damping_0.2_nseq_1.txt" u ($2/1000):9:10 w yerrorl lc 2 lw 2 dt 1 pt 7 notitle \
, "<awk '$1 ==16384 && $3==700 && $4==150' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_seq_directed_ER_summary_T_0.000_lambda_0.000_PD_Lotka_Volterra_final_av0_0.08_dn_0_ninitconds_1_tol_1e-6_maxiter_10000_eps_0.000_damping_0.2_nseq_1.txt" u ($2/1000):9:10 w yerrorl lc 2 lw 2 dt 3 pt 2 notitle \
, "dumb_data.txt" w l lc 2 lw 2 dt 1 title "{/=26 N=16384}" \
, "<awk '$1 ==32768 && $3==700 && $4==0' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_seq_directed_ER_summary_T_0.000_lambda_0.000_PD_Lotka_Volterra_final_av0_0.08_dn_0_ninitconds_1_tol_1e-6_maxiter_10000_eps_0.000_damping_0.2_nseq_1.txt" u ($2/1000):9:10 w yerrorl lc 4 lw 2 dt 1 pt 7 notitle \
, "<awk '$1 ==32768 && $3==700 && $4==150' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_seq_directed_ER_summary_T_0.000_lambda_0.000_PD_Lotka_Volterra_final_av0_0.08_dn_0_ninitconds_1_tol_1e-6_maxiter_10000_eps_0.000_damping_0.2_nseq_1.txt" u ($2/1000):9:10 w yerrorl lc 4 lw 2 dt 3 pt 2 notitle \
, "dumb_data.txt" w l lc 4 lw 2 dt 1 title "{/=26 N=32768}" \
, "<awk '$1 ==65536 && $3==700 && $4==0' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_seq_directed_ER_summary_T_0.000_lambda_0.000_PD_Lotka_Volterra_final_av0_0.08_dn_0_ninitconds_1_tol_1e-6_maxiter_10000_eps_0.000_damping_0.2_nseq_1.txt" u ($2/1000):9:10 w yerrorl lc 6 lw 2 dt 1 pt 7 notitle \
, "<awk '$1 ==65536 && $3==700 && $4==150' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_seq_directed_ER_summary_T_0.000_lambda_0.000_PD_Lotka_Volterra_final_av0_0.08_dn_0_ninitconds_1_tol_1e-6_maxiter_10000_eps_0.000_damping_0.2_nseq_1.txt" u ($2/1000):9:10 w yerrorl lc 6 lw 2 dt 3 pt 2 notitle \
, "dumb_data.txt" w l lc 6 lw 2 dt 1 title "{/=26 N=65536}" \
, 'dumb_data.txt' w yerrorl lc 8 lw 2 dt 1 pt 7 title "{/=26 {/Symbol s} = 0.00}" \
, 'dumb_data.txt' w yerrorl lc 8 lw 2 dt 3 pt 2 title "{/=26 {/Symbol s} = 0.15}"


