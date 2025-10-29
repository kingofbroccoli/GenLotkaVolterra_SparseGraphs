set term postscript enhanced color eps dl 2.5
filenameoutput="IBMF_directed_ER_T0_Lotka_Volterra_N_65536_several_mu.eps"
set output filenameoutput

reset
set multi
set size 1,1
set origin 0,0
set ylabel "{/=26 P_{nc}}" rotate by 90 offset 0,0
set xlabel "{/=26 c}" offset -3,0
set key spacing 2.0 maxrows 10 width 6 at 1.3,0.85
set xra [-0.002:4.0]
set yra [-0.01:1.01]
set tics  font ",24"
set label "{/=28 N = 65536}" at 1.6, 0.3


p "<awk '$1 ==65536 && $3==900 && $4==0' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_T0_directed_ER_PD_Lotka_Volterra_summary_av0_0.08_tol_1e-4_maxiter_10000.txt" u ($2/1000):10:11 w yerrorl lc 1 lw 2 dt 3 pt 7 title "{/=26 {/Symbol m} = 0.90}" \
, "<awk '$1 ==65536 && $3==950 && $4==0' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_T0_directed_ER_PD_Lotka_Volterra_summary_av0_0.08_tol_1e-4_maxiter_10000.txt" u ($2/1000):10:11 w yerrorl lc 2 lw 2 dt 3 pt 4 title "{/=26 {/Symbol m} = 0.95}" \
, "<awk '$1 ==65536 && $3==1050 && $4==0' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_T0_directed_ER_PD_Lotka_Volterra_summary_av0_0.08_tol_1e-4_maxiter_10000.txt" u ($2/1000):10:11 w yerrorl lc 4 lw 2 dt 3 pt 8 title "{/=26 {/Symbol m} = 1.05}" \
, "<awk '$1 ==65536 && $3==1100 && $4==0' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_T0_directed_ER_PD_Lotka_Volterra_summary_av0_0.08_tol_1e-4_maxiter_10000.txt" u ($2/1000):10:11 w yerrorl lc 6 lw 2 dt 3 pt 2 title "{/=26 {/Symbol m} = 1.10}" \
, "<awk '$1 ==65536 && $3==1500 && $4==0' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_T0_directed_ER_PD_Lotka_Volterra_summary_av0_0.08_tol_1e-4_maxiter_10000.txt" u ($2/1000):10:11 w yerrorl lc 7 lw 2 dt 3 pt 10 title "{/=26 {/Symbol m} = 1.50}" \
, 'p_fluc_IBMF_directed_ER_n0_2.txt' w l lc 8 dt 2 title "{/=26 theory}"