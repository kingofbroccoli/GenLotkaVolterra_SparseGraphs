set term postscript enhanced color eps dl 2.5
filenameoutput="Langevin_T0_Lotka_Volterra_crossover.eps"
set output filenameoutput

reset
set size 1,1
set origin 0,0
set ylabel "{/=26 {/Symbol s}}" rotate by 90 offset 0,0
set xlabel "{/=26 {/Symbol m}}" offset 0,0.5
set key spacing 2.0 maxrows 5 width 8 at 0.25,0.16
set xra [-0.01:0.372]
set yra [0.05:0.62]
set tics  font ",24"
set label "{/=26 T=0}" at 0.7, 0.06

set label "{/=26 UNBOUNDED}" at 0.05, 0.56
set label "{/=26 GROWTH}" at 0.065, 0.53

set label "{/=26 SINGLE}" at 0.04, 0.27
set label "{/=26 FIXED}" at 0.047, 0.24
set label "{/=26 POINT}" at 0.046, 0.21


set label "{/=26 MULTIPLE}" at 0.27, 0.45
set label "{/=26 FIXED}" at 0.284, 0.42
set label "{/=26 POINTS}" at 0.279, 0.39



p "../../Data/Simulations/Asymmetric/Lotka-Volterra_transition_mult_epsilon_0.0_Partially_AsymGauss_lambda_1e-06_tol_1e-08_N_1024_c_3.00_T_0.0.txt" u 1:(($2+$3)/2):2:3 w yerror pt 6 lc 6 lw 2 notitle \
, "../../Data/Simulations/Asymmetric/Lotka-Volterra_transition_div_epsilon_0.0_Partially_AsymGauss_lambda_1e-06_tol_1e-08_N_1024_c_3.00_T_0.0.txt" u 1:(($2+$3)/2):2:3 w yerror pt 14 lc 7 lw 2 notitle