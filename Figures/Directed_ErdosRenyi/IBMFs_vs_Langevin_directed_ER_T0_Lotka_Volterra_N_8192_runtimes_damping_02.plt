set term postscript enhanced color eps dl 2.5
filenameoutput="IBMFs_vs_Langevin_directed_ER_T0_Lotka_Volterra_N_8192_runtimes_damping_02.eps"
set output filenameoutput

reset
set multi
set log y
set size 1,1
set origin 0,0
set ylabel "{/=26 runtime(s)}" rotate by 90 offset 0,0
set xlabel "{/=26 c}" offset -3,0
set key spacing 2.0 maxrows 10 width 6 at 3.9,0.1
set xra [-0.002:4.0]
set yra [1e-2:40]
set tics  font ",24"
set label "{/=28 N = 8192}" at 0.6, 0.3


p "<awk '$1 ==8192 && $3==800 && $4==0 && $13-$7>200' ../../Data/Simulations/Directed_ErdosRenyi/Lotka-Volterra_summary_epsilon_0.000_Directed_Trivial_lambda_1e-06_tol_1e-08_T_0.000.txt" u ($2/1000):11:($12/sqrt($13-$7)) w yerror lc 4 lw 3 dt 1 pt 6 title "{/=26 Simulations {/Symbol m} = 0.8}" \
, "<awk '$1 ==8192 && $3==1500 && $4==0 && $13-$7>200' ../../Data/Simulations/Directed_ErdosRenyi/Lotka-Volterra_summary_epsilon_0.000_Directed_Trivial_lambda_1e-06_tol_1e-08_T_0.000.txt" u ($2/1000):11:($12/sqrt($12-$6)) w yerror lc 4 lw 3 dt 1 pt 2 title "{/=26 {/Symbol m} = 1.5}" \
, "<awk '$1 ==8192 && $3==800 && $4==0 && $15-$8>200' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_seq_directed_ER_summary_T_0.000_lambda_0.000_PD_Lotka_Volterra_final_av0_0.08_dn_0_ninitconds_1_tol_1e-6_maxiter_10000_eps_0.000_damping_0.2_nseq_1_Dariah.txt" u ($2/1000):13:($14/sqrt($15-$8)) w yerror lc 6 lw 3 dt 1 pt 4 title "{/=26 IBMF {/Symbol m} = 0.8}"  \
, "<awk '$1 ==8192 && $3==1500 && $4==0 && $15-$8>200' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_seq_directed_ER_summary_T_0.000_lambda_0.000_PD_Lotka_Volterra_final_av0_0.08_dn_0_ninitconds_1_tol_1e-6_maxiter_10000_eps_0.000_damping_0.2_nseq_1_Dariah.txt" u ($2/1000):13:($14/sqrt($15-$8)) w yerror lc 6 lw 3 dt 1 pt 10 title "{/=26 {/Symbol m} = 1.5}"


# reset
# set log y
# set size 0.35,0.35
# set origin 0.1,0.6
# set ylabel "{/=22 Num. Iterations}" rotate by 90 offset 0,0
# set xlabel "{/=22 c}" offset -3,0
# set key spacing 2.0 maxrows 10 width 6 at 1.8,200
# set xra [-0.002:4.0]
# set yra [50:10050]
# set tics  font ",20"
# set xtics 0,1,4


# p "<awk '$1 ==8192 && $3==800 && $4==0 && $13-$7>200' ../../Data/Simulations/Directed_ErdosRenyi/Lotka-Volterra_summary_epsilon_0.000_Directed_Trivial_lambda_1e-06_tol_1e-08_T_0.000.txt" u ($2/1000):5:($6/sqrt($13-$7)) w yerror lc 4 lw 3 dt 1 pt 6 notitle \
# , "<awk '$1 ==8192 && $3==1500 && $4==0 && $13-$7>200' ../../Data/Simulations/Directed_ErdosRenyi/Lotka-Volterra_summary_epsilon_0.000_Directed_Trivial_lambda_1e-06_tol_1e-08_T_0.000.txt" u ($2/1000):5:($6/sqrt($13-$7)) w yerror lc 4 lw 3 dt 1 pt 2 notitle \
# , "<awk '$1 ==8192 && $3==800 && $4==0 && $15-$8>200' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_seq_directed_ER_summary_T_0.000_lambda_0.000_PD_Lotka_Volterra_final_av0_0.08_dn_0_ninitconds_1_tol_1e-6_maxiter_10000_eps_0.000_damping_0.2_nseq_1_Dariah.txt" u ($2/1000):5:($6/sqrt($15-$8)) w yerror lc 6 lw 3 dt 1 pt 4 notitle  \
# , "<awk '$1 ==8192 && $3==1500 && $4==0 && $15-$8>200' ../../Data/IBMF/Directed_ErdosRenyi/IBMF_seq_directed_ER_summary_T_0.000_lambda_0.000_PD_Lotka_Volterra_final_av0_0.08_dn_0_ninitconds_1_tol_1e-6_maxiter_10000_eps_0.000_damping_0.2_nseq_1_Dariah.txt" u ($2/1000):5:($6/sqrt($15-$8)) w yerror lc 6 lw 3 dt 1 pt 10 notitle