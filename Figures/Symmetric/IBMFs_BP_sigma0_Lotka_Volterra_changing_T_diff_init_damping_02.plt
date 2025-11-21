set term postscript enhanced color eps dl 2.5
filenameoutput="IBMFs_BP_sigma0_Lotka_Volterra_changing_T_diff_init_damping_02.eps"
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



p "../../Data/IBMF/Symmetric/IBMF_seq_RRG_PD_Lotka_Volterra_transitions_av0_0.08_lambda_1e-6_tol_1e-6_maxiter_10000_eps_1.000_sigma_0.000_N_1024_c_3_damping_1.0_nseq_10.txt" u (($4+$5)/2):($1):4:5 w xerror lc 6 lw 2 pt 6 title "{/=26 IBMF d=1.0}" \
, "../../Data/IBMF/Symmetric/IBMF_seq_RRG_PD_Lotka_Volterra_transitions_mult_av0_0.5_lambda_1e-6_dn_0.5_ninitconds_10_tol_1e-6_maxiter_10000_eps_1.000_sigma_0.000_N_1024_c_3_damping_0.2_nseq_1.txt" u (($2+$3)/2):($1):2:3 w xerror lc 7 lw 2 pt 2 title "{/=26 IBMF d=0.2}"