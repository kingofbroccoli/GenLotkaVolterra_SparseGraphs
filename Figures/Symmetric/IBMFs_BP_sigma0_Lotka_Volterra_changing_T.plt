set term postscript enhanced color eps dl 2.5
filenameoutput="IBMFs_BP_sigma0_Lotka_Volterra_changing_T.eps"
set output filenameoutput

reset
set multi
set size 1,1
set origin 0,0
set bmargin 3.5
set lmargin 10.5
set ylabel "{/=26 T}" rotate by 90 offset -0.5,0
set xlabel "{/=26 {/Symbol m}}" offset 1,0
set key spacing 2.0 maxrows 10 width 6 at 0.35,0.023
set xra [0.02:0.38]
set yra [-0.0001:0.056]
set tics  font ",24"
set label "{/=26 {/Symbol s}=0}" at 0.04, 0.05



p "../../Data/IBMF/Symmetric/IBMF_seq_RRG_PD_Lotka_Volterra_transitions_av0_0.08_lambda_1e-6_tol_1e-6_maxiter_10000_eps_1.000_sigma_0.000_N_1024_c_3_damping_1.0_nseq_10.txt" u (($4+$5)/2):($1):4:5 w xerror lc 6 lw 2 pt 6 title "{/=26 IBMF}" \
, "../../Data/BP/Symmetric/BP_transition_sigma0_N_128_sequential_changing_order.txt" u ($2-0.0005):($1):($2-0.001):2 w xerror lc 2 lw 2 pt 2 title "{/=26 BP}" \
, "Stars_In_IBMFs_BP_sigma0.txt" u 1:2:3 w p pt 59 ps 1.5 lc variable lw 2 notitle