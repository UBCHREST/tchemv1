#
#   Ignition delay time
#
unset label
# set nokey
set size 0.75,1
set style line 1 lt  1 lw 3
set style line 2 lt  2 lw 3
set style line 3 lt  3 lw 3
set style line 4 lt  4 lw 3
set style line 5 lt -1 lw 3
set border 15 lw 2
set xtics font "Times,18"
set ytics font "Times,18"
set xlabel "1000/T" font "Times,24"
set ylabel "Ignition delay time [s]" font "Times,24"
set log y
set term postscript eps color enhanced
set output "gri3_igndel.eps"
set lmargin 12 
set rmargin 3 
set tmargin 2 
set bmargin 4 
p 'igndel_2_t.dat'  u (1000/$1):2 w lp lt  1 lw 3 pt 1 ps 1 title "2 atm", \
  'igndel_10_t.dat' u (1000/$1):2 w lp lt  2 lw 3 pt 1 ps 1 title "10 atm", \
  'igndel_50_t.dat' u (1000/$1):2 w lp lt  3 lw 3 pt 1 ps 1 title "50 atm"


