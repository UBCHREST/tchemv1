#
#   Ignition delay time
#
unset label
set key right bottom
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
set output "iOct_id2.eps"
set lmargin 12 
set rmargin 3 
set tmargin 2 
set bmargin 4 
p 'igndel_13_t.dat' u (1000/$1):2 w lp lt  1 lw 3 pt 1 ps 1 title "13 atm", \
  'igndel_16_t.dat' u (1000/$1):2 w lp lt  2 lw 3 pt 2 ps 1 title "16 atm", \
  'igndel_34_t.dat' u (1000/$1):2 w lp lt  3 lw 3 pt 3 ps 1 title "34 atm", \
  'igndel_40_t.dat' u (1000/$1):2 w lp lt  4 lw 3 pt 4 ps 1 title "40 atm", \
  'igndel_45_t.dat' u (1000/$1):2 w lp lt -1 lw 3 pt 5 ps 1 title "45 atm"


