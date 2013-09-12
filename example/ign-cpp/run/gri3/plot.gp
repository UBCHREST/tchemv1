#
#   Major species: O2,CH4,H2O,CO2,CO
#
unset label
set nokey
set size 0.75,1
set style line 1 lt  1 lw 8
set style line 2 lt  2 lw 8
set style line 3 lt  3 lw 8
set style line 4 lt  4 lw 8
set style line 5 lt -1 lw 4
set border 15 lw 2
set xtics font "Times,18"
set ytics font "Times,18"
set xlabel "time [s]" font "Times,24"
set ylabel "Mass fraction" font "Times,24"
set label "O_2"   at graph 0.20,0.88 font "Times,18"
set label "CH_4"  at graph 0.20,0.22 font "Times,18"
set label "H_2O"  at graph 0.80,0.52 font "Times,18"
set label "CO_2"  at graph 0.80,0.42 font "Times,18"
set label "CO"    at graph 0.80,0.16 font "Times,18"
set term postscript eps color enhanced
set output "specMaj.eps"
set lmargin 10 
set rmargin 3 
set tmargin 2 
set bmargin 4 
p [1.099:1.102] [0:2.3e-1] 'ys.out' u 1:6  w l ls 1,\
                       'ys.out' u 1:16 w l ls 2,\
                       'ys.out' u 1:18 w l ls 3,\
                       'ys.out' u 1:8  w l ls 4,\
                       'ys.out' u 1:17 w l ls 5


#
#  Minor species CH3, HCO, CH2O, C2H2
#
reset
unset label
set nokey
set size 0.75,1
set style line 1 lt  1 lw 8
set style line 2 lt  2 lw 8
set style line 3 lt  3 lw 8
set style line 4 lt  4 lw 8
set style line 5 lt -1 lw 4
set border 15 lw 2
set xtics font "Times,18"
set ytics font "Times,18"
set xlabel "time [s]"      font "Times,24"
set ylabel "Mass fraction" font "Times,24"
set label "CH_2(*20)"  at graph 0.20,0.07 font "Times,18"
set label "CH(*200)"   at graph 0.52,0.85 font "Times,18"
set label "HCO(*10)"   at graph 0.20,0.17 font "Times,18"
set label "C_2H_2"     at graph 0.20,0.30 font "Times,18"
set lmargin 12 
set rmargin 4 
set tmargin 2 
set bmargin 4 
#set arrow 1 from graph 0.82,0.70 to graph 0.64,0.50 lt -1
#set arrow 2 from graph 0.82,0.50 to graph 0.68,0.42 lt -1
set term postscript eps color enhanced
set xtics 1.1008,0.00003
set output "specMin.eps"
p [1.1008:1.10088] [0:1.1e-3] \
 'ys.out' u 1:($13*20)  w l ls 1, \
 'ys.out' u 1:($19*10)  w l ls 2, \
 'ys.out' u 1:($12*200) w l ls 3, \
 'ys.out' u 1:($25*1)   w l ls 5   


