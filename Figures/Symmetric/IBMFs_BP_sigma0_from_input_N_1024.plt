set term postscript enhanced color eps dl 2.5
filenameoutput="IBMFs_BP_sigma0_from_input_N_1024.eps"
set output filenameoutput

reset
set multi
set bmargin 4
set size 1,1
set origin 0,0
set ylabel "{/=26 <n>}" rotate by 90 offset 0,0
set xlabel "{/=26 {/Symbol m}}" offset 0,0
set key spacing 2.0 maxrows 3 width 6 at 0.3,0.35
set xra [0:0.356]
set yra [0.16:1.1]
set tics  font ",24"

p "<awk '$1==1' ../../Data/IBMF/Symmetric/IBMF_seq_from_input_RRG_PD_Lotka_Volterra_summary_av0_0.08_lambda_1e-6_tol_1e-6_maxiter_10000_eps_1.0_sigma_0.000_N_1024_c_3_damping_1.0_nseq_10000.txt" u ($2<349?$2/1000:1/0):8 w l lc 1 lw 3 dt 1 notitle \
, "<awk '$1==10' ../../Data/IBMF/Symmetric/IBMF_seq_from_input_RRG_PD_Lotka_Volterra_summary_av0_0.08_lambda_1e-6_tol_1e-6_maxiter_10000_eps_1.0_sigma_0.000_N_1024_c_3_damping_1.0_nseq_10000.txt" u ($2<243?$2/1000:1/0):8 w l lc 2 lw 3 dt 1 notitle \
, "<awk '$1==20' ../../Data/IBMF/Symmetric/IBMF_seq_from_input_RRG_PD_Lotka_Volterra_summary_av0_0.08_lambda_1e-6_tol_1e-6_maxiter_10000_eps_1.0_sigma_0.000_N_1024_c_3_damping_1.0_nseq_10000.txt" u ($2<141?$2/1000:1/0):8 w l lc 3 lw 3 dt 1 notitle \
, "<awk '$1==25' ../../Data/IBMF/Symmetric/IBMF_seq_from_input_RRG_PD_Lotka_Volterra_summary_av0_0.08_lambda_1e-6_tol_1e-6_maxiter_10000_eps_1.0_sigma_0.000_N_1024_c_3_damping_1.0_nseq_10000.txt" u ($2<101?$2/1000:1/0):8 w l lc 4 lw 3 dt 1 notitle \
, "<awk '$1==28' ../../Data/IBMF/Symmetric/IBMF_seq_from_input_RRG_PD_Lotka_Volterra_summary_av0_0.08_lambda_1e-6_tol_1e-6_maxiter_10000_eps_1.0_sigma_0.000_N_1024_c_3_damping_1.0_nseq_10000.txt" u ($2<85?$2/1000:1/0):8 w l lc 5 lw 3 dt 1 notitle \
, "<awk '$1==31' ../../Data/IBMF/Symmetric/IBMF_seq_from_input_RRG_PD_Lotka_Volterra_summary_av0_0.08_lambda_1e-6_tol_1e-6_maxiter_10000_eps_1.0_sigma_0.000_N_1024_c_3_damping_1.0_nseq_10000.txt" u ($2<72?$2/1000:1/0):8 w l lc 7 lw 3 dt 1 notitle \
, "<awk '$1==34' ../../Data/IBMF/Symmetric/IBMF_seq_from_input_RRG_PD_Lotka_Volterra_summary_av0_0.08_lambda_1e-6_tol_1e-6_maxiter_10000_eps_1.0_sigma_0.000_N_1024_c_3_damping_1.0_nseq_10000.txt" u ($2<65?$2/1000:1/0):8 w l lc 6 lw 3 dt 1 notitle \
, "<awk '$1==1' ../../Data/IBMF/Symmetric/IBMF_seq_from_input_RRG_PD_Lotka_Volterra_summary_av0_0.08_lambda_1e-6_tol_1e-6_maxiter_10000_eps_1.0_sigma_0.000_N_1024_c_3_damping_1.0_nseq_10000.txt" u ($2>349?$2/1000:1/0):8:9 w yerrorl lc 1 lw 3 dt 1 pt 2 ps 1.2 title "{/=26 T=0.001}" \
, "<awk '$1==10' ../../Data/IBMF/Symmetric/IBMF_seq_from_input_RRG_PD_Lotka_Volterra_summary_av0_0.08_lambda_1e-6_tol_1e-6_maxiter_10000_eps_1.0_sigma_0.000_N_1024_c_3_damping_1.0_nseq_10000.txt" u ($2>243?$2/1000:1/0):8:9 w yerrorl lc 2 lw 3 dt 1 pt 4 ps 1.2 title "{/=26 T=0.010}" \
, "<awk '$1==20' ../../Data/IBMF/Symmetric/IBMF_seq_from_input_RRG_PD_Lotka_Volterra_summary_av0_0.08_lambda_1e-6_tol_1e-6_maxiter_10000_eps_1.0_sigma_0.000_N_1024_c_3_damping_1.0_nseq_10000.txt" u ($2>141?$2/1000:1/0):8:9 w yerrorl lc 3 lw 3 dt 1 pt 6 ps 1.2 title "{/=26 T=0.020}" \
, "<awk '$1==25' ../../Data/IBMF/Symmetric/IBMF_seq_from_input_RRG_PD_Lotka_Volterra_summary_av0_0.08_lambda_1e-6_tol_1e-6_maxiter_10000_eps_1.0_sigma_0.000_N_1024_c_3_damping_1.0_nseq_10000.txt" u ($2>101?$2/1000:1/0):8:9 w yerrorl lc 4 lw 3 dt 1 pt 8 ps 1.2 title "{/=26 T=0.025}" \
, "<awk '$1==28' ../../Data/IBMF/Symmetric/IBMF_seq_from_input_RRG_PD_Lotka_Volterra_summary_av0_0.08_lambda_1e-6_tol_1e-6_maxiter_10000_eps_1.0_sigma_0.000_N_1024_c_3_damping_1.0_nseq_10000.txt" u ($2>85?$2/1000:1/0):8:9 w yerrorl lc 5 lw 3 dt 1 pt 10 ps 1.2 title "{/=26 T=0.028}" \
, "<awk '$1==31' ../../Data/IBMF/Symmetric/IBMF_seq_from_input_RRG_PD_Lotka_Volterra_summary_av0_0.08_lambda_1e-6_tol_1e-6_maxiter_10000_eps_1.0_sigma_0.000_N_1024_c_3_damping_1.0_nseq_10000.txt" u ($2>72?$2/1000:1/0):8:9 w yerrorl lc 7 lw 3 dt 1 pt 12 ps 1.2 title "{/=26 T=0.031}" \
, "<awk '$1==34' ../../Data/IBMF/Symmetric/IBMF_seq_from_input_RRG_PD_Lotka_Volterra_summary_av0_0.08_lambda_1e-6_tol_1e-6_maxiter_10000_eps_1.0_sigma_0.000_N_1024_c_3_damping_1.0_nseq_10000.txt" u ($2>65?$2/1000:1/0):8:9 w yerrorl lc 6 lw 3 dt 1 pt 14 ps 1.2 title "{/=26 T=0.034}" \
, "<awk '$1==1' ../../Data/BP/Symmetric/BP_from_input_RRG_PD_Lotka_Volterra_summary_N_1024_damping_1.0.txt" u ($2/1000):7 w l lc 8 lw 4 dt 3 title "{/=26 BP}" \
, "<awk '$1==10' ../../Data/BP/Symmetric/BP_from_input_RRG_PD_Lotka_Volterra_summary_N_1024_damping_1.0.txt" u ($2/1000):7 w l lc 8 lw 4 dt 3 notitle \
, "<awk '$1==20' ../../Data/BP/Symmetric/BP_from_input_RRG_PD_Lotka_Volterra_summary_N_1024_damping_1.0.txt" u ($2/1000):7 w l lc 8 lw 4 dt 3 notitle \
, "<awk '$1==25' ../../Data/BP/Symmetric/BP_from_input_RRG_PD_Lotka_Volterra_summary_N_1024_damping_1.0.txt" u ($2/1000):7 w l lc 8 lw 4 dt 3 notitle \
, "<awk '$1==28' ../../Data/BP/Symmetric/BP_from_input_RRG_PD_Lotka_Volterra_summary_N_1024_damping_1.0.txt" u ($2/1000):7 w l lc 8 lw 4 dt 3 notitle \
, "<awk '$1==31' ../../Data/BP/Symmetric/BP_from_input_RRG_PD_Lotka_Volterra_summary_N_1024_damping_1.0.txt" u ($2/1000):7 w l lc 8 lw 4 dt 3 notitle \
, "<awk '$1==34' ../../Data/BP/Symmetric/BP_from_input_RRG_PD_Lotka_Volterra_summary_N_1024_damping_1.0.txt" u ($2/1000):7 w l lc 8 lw 4 dt 3 notitle

reset

set size 0.45,0.45
set origin 0.5,0.52
set ylabel "{/=22 {/Symbol g} (X 10^{4})}" rotate by 90 offset -1,0
set xlabel "{/=22 {/Symbol m}}" offset 0,0
set key spacing 2.0 maxrows 3 width 6 at 0.3,0.35
set xra [0:0.356]
set yra [-6:4]
set tics  font ",20"
set ytics -5,4,6
set xtics 0,0.1,0.3

p "<awk '$1==1' ../../Data/BP/Symmetric/BP_from_input_RRG_PD_Lotka_Volterra_summary_N_1024_damping_1.0.txt" u ($2<349?$2/1000:1/0):($11 * 10000) w l lc 1 lw 4 dt 3 notitle \
, "<awk '$1==10' ../../Data/BP/Symmetric/BP_from_input_RRG_PD_Lotka_Volterra_summary_N_1024_damping_1.0.txt" u ($2<243?$2/1000:1/0):($11 * 10000) w l lc 2 lw 4 dt 3 notitle \
, "<awk '$1==20' ../../Data/BP/Symmetric/BP_from_input_RRG_PD_Lotka_Volterra_summary_N_1024_damping_1.0.txt" u ($2<141?$2/1000:1/0):($11 * 10000) w l lc 3 lw 4 dt 3 notitle \
, "<awk '$1==25' ../../Data/BP/Symmetric/BP_from_input_RRG_PD_Lotka_Volterra_summary_N_1024_damping_1.0.txt" u ($2<101?$2/1000:1/0):($11 * 10000) w l lc 4 lw 4 dt 3 notitle \
, "<awk '$1==28' ../../Data/BP/Symmetric/BP_from_input_RRG_PD_Lotka_Volterra_summary_N_1024_damping_1.0.txt" u ($2<85?$2/1000:1/0):($11 * 10000) w l lc 5 lw 4 dt 3 notitle \
, "<awk '$1==31' ../../Data/BP/Symmetric/BP_from_input_RRG_PD_Lotka_Volterra_summary_N_1024_damping_1.0.txt" u ($2<72?$2/1000:1/0):($11 * 10000) w l lc 7 lw 4 dt 3 notitle \
, "<awk '$1==34' ../../Data/BP/Symmetric/BP_from_input_RRG_PD_Lotka_Volterra_summary_N_1024_damping_1.0.txt" u ($2<65?$2/1000:1/0):($11 * 10000) w l lc 6 lw 4 dt 3 notitle \
, "<awk '$1==1' ../../Data/BP/Symmetric/BP_from_input_RRG_PD_Lotka_Volterra_summary_N_1024_damping_1.0.txt" u ($2>349?$2/1000:1/0):($11 * 10000) w l lc 8 lw 4 dt 1 notitle \
, "<awk '$1==10' ../../Data/BP/Symmetric/BP_from_input_RRG_PD_Lotka_Volterra_summary_N_1024_damping_1.0.txt" u ($2>243?$2/1000:1/0):($11 * 10000) w l lc 8 lw 3 dt 1 notitle \
, "<awk '$1==20' ../../Data/BP/Symmetric/BP_from_input_RRG_PD_Lotka_Volterra_summary_N_1024_damping_1.0.txt" u ($2>141?$2/1000:1/0):($11 * 10000) w l lc 8 lw 3 dt 1 notitle \
, "<awk '$1==25' ../../Data/BP/Symmetric/BP_from_input_RRG_PD_Lotka_Volterra_summary_N_1024_damping_1.0.txt" u ($2>101?$2/1000:1/0):($11 * 10000) w l lc 8 lw 3 dt 1 notitle \
, "<awk '$1==28' ../../Data/BP/Symmetric/BP_from_input_RRG_PD_Lotka_Volterra_summary_N_1024_damping_1.0.txt" u ($2>85?$2/1000:1/0):($11 * 10000) w l lc 8 lw 3 dt 1 notitle \
, "<awk '$1==31' ../../Data/BP/Symmetric/BP_from_input_RRG_PD_Lotka_Volterra_summary_N_1024_damping_1.0.txt" u ($2>72?$2/1000:1/0):($11 * 10000) w l lc 8 lw 3 dt 1 notitle \
, "<awk '$1==34' ../../Data/BP/Symmetric/BP_from_input_RRG_PD_Lotka_Volterra_summary_N_1024_damping_1.0.txt" u ($2>65?$2/1000:1/0):($11 * 10000) w l lc 8 lw 3 dt 1 notitle