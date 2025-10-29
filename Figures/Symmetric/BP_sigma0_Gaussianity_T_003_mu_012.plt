set term postscript enhanced color eps dl 2.5
set encoding utf8
filenameoutput="BP_sigma0_Gaussianity_T_003_mu_012.eps"
set output filenameoutput

reset
set multi
set size 1,1
set origin 0,0
set log y
set bmargin 3.5
set lmargin 12.5
set ylabel "{/=26 ~P{.55\U+007E}_{BP}(n)}" rotate by 90 offset -1,2
set xlabel "{/=26 n}" offset 1,0
set key spacing 2.0 maxrows 10 width 6 at 1.4,0.03
set yrange [1e-5:2.5]
set xrange [0.00:1.87]
set tics  font ",24"
set format y "%g"

set label "" at -0.2, 1.6 point pt 59 lc 7 ps 4

p "../../Data/BP/Symmetric/BP_cont_LV_sigma0_true_marginals_gaussian_part_N_128_T_0.030_mu_0.120.txt" u 1:2 w p pt 10 lw 2 lc 7 ps 1.2 title "{/=26 BP  T=0.03  {/Symbol m}=0.12}" \
, "../../Data/BP/Symmetric/BP_cont_LV_sigma0_true_marginals_fit_trunc_gauss_N_128_T_0.030_mu_0.120.txt" u 1:2 w l lc 8 lw 2 dt 1 title "{/=26 Fitted Trunc. Gauss.}"


reset
set size 0.4,0.4
set origin 0.32,0.1
set yrange [-0.6:1]
set xrange [0.06:1.76]
set bmargin 3.5
set lmargin 10.5
set ylabel "{/=20 {/Symbol D}_{BP-G}}" rotate by 90 offset -0.4,0
set xlabel "{/=20 n}" offset 1.8,0.8
set key spacing 1.5 maxrows 5 width 3 at 3.2,0.9
set tics  font ",18"
set xtics 0.4,0.4,1.6
set ytics -0.4,0.4,1.0


p "../../Data/BP/Symmetric/BP_cont_LV_sigma0_diff_true_marginals_fit_trunc_gauss_N_128_T_0.030_mu_0.120.txt" ev 5:5 w p pt 10 lw 2 lc 7 ps 0.5 notitle