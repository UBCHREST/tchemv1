#
#   Ignition delay time
#
unset label
# set nokey
set size 0.75,1
set style line 1 lt  1 lw 8
set style line 2 lt  2 lw 8
set style line 3 lt  3 lw 8
set style line 4 lt  4 lw 8
set style line 5 lt -1 lw 8
set border 15 lw 2
set xtics font "Times,18"
set ytics font "Times,18"
set xlabel "f" font "Symbol,24"
set ylabel "Ignition delay time [s]" font "Times,24"
set log y
#set label "O_2"   at graph 0.20,0.88 font "Times,18"
#set label "CH_4"  at graph 0.20,0.22 font "Times,18"
set term postscript eps color enhanced
set output "igndel.eps"
set lmargin 10 
set rmargin 3 
set tmargin 2 
set bmargin 4 
p [0:3.2] 'igndel_1atm.dat'  u 1:2 title "1 atm"  w l ls 1, \
          'igndel_10atm.dat' u 1:2 title "10 atm" w l ls 2, \
          'igndel_20atm.dat' u 1:2 title "20 atm" w l ls 3

