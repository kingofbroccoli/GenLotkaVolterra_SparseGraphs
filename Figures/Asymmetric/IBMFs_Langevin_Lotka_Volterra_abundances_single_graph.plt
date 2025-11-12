set term postscript enhanced color eps dl 2.5
filenameoutput="IBMFs_Langevin_Lotka_Volterra_abundances_single_graph.eps"
set output filenameoutput

reset
set multi
set size 1,1
set origin 0,0
set bmargin 3.2
set ylabel "{/=26 n^{IBMF}}" rotate by 90 offset 1,0
set xlabel "{/=26 n^{SIM}}" offset 6,0.5
set key spacing 2.0 maxrows 5 width 8 at 0.29,0.13
set xra [0:2]
set yra [0:2]
set tics  font ",24"


id1=303
id2=194
id3=396
id4=149

id5=78
id6=699
id7=398
id8=99


p "IBMF_vs_Langevin_single_abundances.txt" u 1:4 w p pt 4 lc 6 lw 2 notitle \
, sprintf("< awk '$3 == %d' IBMF_vs_Langevin_single_abundances.txt", id1) u 1:4 w p pt 5 lc 1 lw 3 notitle \
, sprintf("< awk '$3 == %d' IBMF_vs_Langevin_single_abundances.txt", id2) u 1:4 w p pt 5 lc 2 lw 3 notitle \
, sprintf("< awk '$3 == %d' IBMF_vs_Langevin_single_abundances.txt", id3) u 1:4 w p pt 5 lc 7 lw 3 notitle \
, sprintf("< awk '$3 == %d' IBMF_vs_Langevin_single_abundances.txt", id4) u 1:4 w p pt 5 lc 4 lw 3 notitle \
, sprintf("< awk '$3 == %d' IBMF_vs_Langevin_single_abundances.txt", id5) u 1:4 w p pt 5 lc rgb "#4a0c6b" lw 3 notitle \
, sprintf("< awk '$3 == %d' IBMF_vs_Langevin_single_abundances.txt", id6) u 1:4 w p pt 5 lc rgb "#6ece58" lw 3 notitle \
, sprintf("< awk '$3 == %d' IBMF_vs_Langevin_single_abundances.txt", id7) u 1:4 w p pt 5 lc rgb "#f7d13d" lw 3 notitle \
, sprintf("< awk '$3 == %d' IBMF_vs_Langevin_single_abundances.txt", id8) u 1:4 w p pt 5 lc rgb "#bd3786" lw 3 notitle \
, sprintf("< awk '$3 == %d' IBMF_vs_Langevin_single_abundances.txt", id1) u 1:4 w p pt 4 lc 8 lw 3 notitle \
, sprintf("< awk '$3 == %d' IBMF_vs_Langevin_single_abundances.txt", id2) u 1:4 w p pt 4 lc 8 lw 3 notitle \
, sprintf("< awk '$3 == %d' IBMF_vs_Langevin_single_abundances.txt", id3) u 1:4 w p pt 4 lc 8 lw 3 notitle \
, sprintf("< awk '$3 == %d' IBMF_vs_Langevin_single_abundances.txt", id4) u 1:4 w p pt 4 lc 8 lw 3 notitle \
, sprintf("< awk '$3 == %d' IBMF_vs_Langevin_single_abundances.txt", id5) u 1:4 w p pt 4 lc 8 lw 3 notitle \
, sprintf("< awk '$3 == %d' IBMF_vs_Langevin_single_abundances.txt", id6) u 1:4 w p pt 4 lc 8 lw 3 notitle \
, sprintf("< awk '$3 == %d' IBMF_vs_Langevin_single_abundances.txt", id7) u 1:4 w p pt 4 lc 8 lw 3 notitle \
, sprintf("< awk '$3 == %d' IBMF_vs_Langevin_single_abundances.txt", id8) u 1:4 w p pt 4 lc 8 lw 3 notitle \
, x w l lc 8 lw 2 dt 2 notitle


reset
set log y
set size 0.45,0.45
set origin 0.07,0.53
set bmargin 3.2
set ylabel "{/=22 n_{i}^{SIM}}" rotate by 90 offset 4,0
set xlabel "{/=22 t}" offset 0,0
set key spacing 1.5 maxrows 5 width 5 at 780,1.13
set xra [0:500]
set yra [1e-5:5]
set tics  font ",20"
set xtics 0, 250, 500
set ytics add ("" 1e-5, "1e-4" 1e-4, "" 1e-3, "0.01" 1e-2, "" 1e-1, "1" 1)
set mytics default
set label "{/=22 Small n_i}" at 150,1.2



val1 = system(sprintf("awk '($3==%d) {print $4}' IBMF_vs_Langevin_single_abundances.txt", id1))
val1 = real(val1)
val2 = system(sprintf("awk '($3==%d) {print $4}' IBMF_vs_Langevin_single_abundances.txt", id2))
val2 = real(val2)
val3 = system(sprintf("awk '($3==%d) {print $4}' IBMF_vs_Langevin_single_abundances.txt", id3))
val3 = real(val3)
val4 = system(sprintf("awk '($3==%d) {print $4}' IBMF_vs_Langevin_single_abundances.txt", id4))
val4 = real(val4)

col1 = id1 + 2
col2 = id2 + 2
col3 = id3 + 2
col4 = id4 + 2

p "../../Langevin/Results/One_graph/Lotka-Volterra_Extraction_1_Measure_2_mu_0.0_sigma_0.150_T_0.015.txt" u 1:col1 w l lc 1 lw 2 notitle\
, "../../Langevin/Results/One_graph/Lotka-Volterra_Extraction_1_Measure_2_mu_0.0_sigma_0.150_T_0.015.txt" u 1:col2 w l lc 2 lw 2 notitle \
, "../../Langevin/Results/One_graph/Lotka-Volterra_Extraction_1_Measure_2_mu_0.0_sigma_0.150_T_0.015.txt" u 1:col3 w l lc 7 lw 2 notitle \
, "../../Langevin/Results/One_graph/Lotka-Volterra_Extraction_1_Measure_2_mu_0.0_sigma_0.150_T_0.015.txt" u 1:col4 w l lc 4 lw 2 notitle \
, val1 w l lc 1 lw 3 dt 3 notitle \
, val2 w l lc 2 lw 3 dt 3 notitle \
, val3 w l lc 7 lw 3 dt 3 notitle \
, val4 w l lc 4 lw 3 dt 3 notitle



reset
set size 0.45,0.45
set origin 0.5,0.08
set bmargin 3.2
set ylabel "{/=22 n_{i}^{SIM}}" rotate by 90 offset 0,0
set xlabel "{/=22 t}" offset 4,1
set key spacing 1.5 maxrows 5 width 5 at 780,1.13
set xra [0:500]
set yra [0.4:2]
set tics  font ",20"
set xtics 0, 250, 500
set ytics 0.4,0.4,2
set label "{/=22 Large n_i}" at 160,0.55


val5 = system(sprintf("awk '($3==%d) {print $4}' IBMF_vs_Langevin_single_abundances.txt", id5))
val5 = real(val5)
val6 = system(sprintf("awk '($3==%d) {print $4}' IBMF_vs_Langevin_single_abundances.txt", id6))
val6 = real(val6)
val7 = system(sprintf("awk '($3==%d) {print $4}' IBMF_vs_Langevin_single_abundances.txt", id7))
val7 = real(val7)
val8 = system(sprintf("awk '($3==%d) {print $4}' IBMF_vs_Langevin_single_abundances.txt", id8))
val8 = real(val8)

col5 = id5 + 2
col6 = id6 + 2
col7 = id7 + 2
col8 = id8 + 2

p "../../Langevin/Results/One_graph/Lotka-Volterra_Extraction_1_Measure_2_mu_0.0_sigma_0.150_T_0.015.txt" u 1:col5 w l lc rgb "#4a0c6b" lw 2 notitle\
, "../../Langevin/Results/One_graph/Lotka-Volterra_Extraction_1_Measure_2_mu_0.0_sigma_0.150_T_0.015.txt" u 1:col6 w l lc rgb "#6ece58" lw 2 notitle \
, "../../Langevin/Results/One_graph/Lotka-Volterra_Extraction_1_Measure_2_mu_0.0_sigma_0.150_T_0.015.txt" u 1:col7 w l lc rgb "#f7d13d" lw 2 notitle \
, "../../Langevin/Results/One_graph/Lotka-Volterra_Extraction_1_Measure_2_mu_0.0_sigma_0.150_T_0.015.txt" u 1:col8 w l lc rgb "#bd3786" lw 2 notitle \
, val5 w l lc 8 lw 3 dt 3 notitle \
, val6 w l lc 8 lw 3 dt 3 notitle \
, val7 w l lc 8 lw 3 dt 3 notitle \
, val8 w l lc 8 lw 3 dt 3 notitle