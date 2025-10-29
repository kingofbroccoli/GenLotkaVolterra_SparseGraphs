set term postscript enhanced color eps dl 2.5
filenameoutput="IBMF_directed_ER_T0_Lotka_Volterra_several_N.eps"
set output filenameoutput

reset
set multi
set size 1,1
set origin 0,0
set ylabel "{/=26 P_{nc}}" rotate by 90 offset 0,0
set xlabel "{/=26 c}" offset -3,0
set key spacing 2.0 maxrows 10 width 6 at 3.9,0.47
set xra [-0.002:4.0]
set yra [-0.01:1.01]
set tics  font ",24"


p "<awk '$1 ==8192 && $3==900 && $4==0' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_T0_directed_ER_PD_Lotka_Volterra_summary_av0_0.08_tol_1e-4_maxiter_10000.txt" u ($2/1000):10:11 w yerrorl lc 7 lw 2 dt 3 pt 7 notitle \
, "<awk '$1 ==8192 && $3==1100 && $4==0' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_T0_directed_ER_PD_Lotka_Volterra_summary_av0_0.08_tol_1e-4_maxiter_10000.txt" u ($2/1000):10:11 w yerrorl lc 7 lw 2 dt 3 pt 2 notitle \
, "dumb_data.txt" w l lc 7 lw 2 dt 3 title "{/=26 N=8192}" \
, "<awk '$1 ==16384 && $3==900 && $4==0' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_T0_directed_ER_PD_Lotka_Volterra_summary_av0_0.08_tol_1e-4_maxiter_10000.txt" u ($2/1000):10:11 w yerrorl lc 2 lw 2 dt 3 pt 7 notitle \
, "<awk '$1 ==16384 && $3==1100 && $4==0' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_T0_directed_ER_PD_Lotka_Volterra_summary_av0_0.08_tol_1e-4_maxiter_10000.txt" u ($2/1000):10:11 w yerrorl lc 2 lw 2 dt 3 pt 2 notitle \
, "dumb_data.txt" w l lc 2 lw 2 dt 3 title "{/=26 N=16384}" \
, "<awk '$1 ==32768 && $3==900 && $4==0' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_T0_directed_ER_PD_Lotka_Volterra_summary_av0_0.08_tol_1e-4_maxiter_10000.txt" u ($2/1000):10:11 w yerrorl lc 4 lw 2 dt 3 pt 7 notitle \
, "<awk '$1 ==32768 && $3==1100 && $4==0' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_T0_directed_ER_PD_Lotka_Volterra_summary_av0_0.08_tol_1e-4_maxiter_10000.txt" u ($2/1000):10:11 w yerrorl lc 4 lw 2 dt 3 pt 2 notitle \
, "dumb_data.txt" w l lc 4 lw 2 dt 3 title "{/=26 N=32768}" \
, "<awk '$1 ==65536 && $3==900 && $4==0' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_T0_directed_ER_PD_Lotka_Volterra_summary_av0_0.08_tol_1e-4_maxiter_10000.txt" u ($2/1000):10:11 w yerrorl lc 6 lw 2 dt 3 pt 7 notitle \
, "<awk '$1 ==65536 && $3==1100 && $4==0' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_T0_directed_ER_PD_Lotka_Volterra_summary_av0_0.08_tol_1e-4_maxiter_10000.txt" u ($2/1000):10:11 w yerrorl lc 6 lw 2 dt 3 pt 2 notitle \
, "dumb_data.txt" w l lc 6 lw 2 dt 3 title "{/=26 N=65536}" \
, 'dumb_data.txt' w yerrorl lc 8 lw 2 dt 3 pt 7 title "{/=26 {/Symbol m} = 0.9}" \
, 'dumb_data.txt' w yerrorl lc 8 lw 2 dt 3 pt 2 title "{/=26 {/Symbol m} = 1.1}" \
, 'p_fluc_IBMF_directed_ER_n0_2.txt' w l lc 8 dt 2 title "{/=26 theory}"



reset
set size 0.43,0.43
set origin 0.07,0.5
set ylabel "{/=22 P_{conv}}" rotate by 90 offset 4,0
set xlabel "{/=22 c}" offset 1.6,0.2
set key spacing 2.0 maxrows 10 width 6 at 3.45,0.4
set xra [2.65:2.8]
set yra [0.95:1.01]
set tics  font ",20"
set ytics (0.95, 1)
set xtics (2.68, "e" exp(1), 2.76)
set arrow nohead from exp(1), 0.95 to exp(1),1.01 dt 2 lc 8


p "<awk '$1 ==8192 && $3==1100 && $4==0' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_T0_directed_ER_PD_Lotka_Volterra_summary_av0_0.08_tol_1e-4_maxiter_10000.txt" u ($2/1000):10:11 w yerror lc 7 lw 2 dt 3 pt 2 notitle \
, "<awk '$1 ==16384 && $3==1100 && $4==0' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_T0_directed_ER_PD_Lotka_Volterra_summary_av0_0.08_tol_1e-4_maxiter_10000.txt" u ($2/1000):10:11 w yerror lc 2 lw 2 dt 3 pt 2 notitle \
, "<awk '$1 ==32768 && $3==1100 && $4==0' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_T0_directed_ER_PD_Lotka_Volterra_summary_av0_0.08_tol_1e-4_maxiter_10000.txt" u ($2/1000):10:11 w yerror lc 4 lw 2 dt 3 pt 2 notitle \
, "<awk '$1 ==65536 && $3==1100 && $4==0' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_T0_directed_ER_PD_Lotka_Volterra_summary_av0_0.08_tol_1e-4_maxiter_10000.txt" u ($2/1000):10:11 w yerror lc 6 lw 2 dt 3 pt 2 notitle