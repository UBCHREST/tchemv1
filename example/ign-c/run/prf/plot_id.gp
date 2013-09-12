#
#   Ignition delay time
#
unset label
# set nokey
set key left top
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
set output "prf_id1.eps"
set lmargin 12 
set rmargin 3 
set tmargin 2 
set bmargin 4 
p 'igndel_p15atm_on100_t.dat' u (1000/$1):2 w lp lt 1 lw 3 pt 1 ps 1 pointinterval 4 title "15 atm, ON 100", \
  'igndel_p15atm_on80_t.dat'  u (1000/$1):2 w lp lt 1 lw 3           pointinterval 4 title "15 atm, ON   80", \
  'igndel_p15atm_on60_t.dat'  u (1000/$1):2 w lp lt 1 lw 3           pointinterval 4 title "15 atm, ON   60", \
  'igndel_p45atm_on100_t.dat' u (1000/$1):2 w lp lt 3 lw 3 pt 4 ps 1 pointinterval 4 title "45 atm, ON 100", \
  'igndel_p45atm_on80_t.dat'  u (1000/$1):2 w lp lt 3 lw 3           pointinterval 4 title "45 atm, ON   80", \
  'igndel_p45atm_on60_t.dat'  u (1000/$1):2 w lp lt 3 lw 3           pointinterval 4 title "45 atm, ON   60"

#
#   Ignition delay time
#
unset label
# set nokey
set key left top
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
set output "prf_id2.eps"
set lmargin 12 
set rmargin 3 
set tmargin 2 
set bmargin 4 
p 'igndel_p15atm_on50_t.dat' u (1000/$1):2 w l lt 1 lw 3 title "15 atm", \
  'igndel_p30atm_on50_t.dat' u (1000/$1):2 w l lt 2 lw 3 title "30 atm", \
  'igndel_p45atm_on50_t.dat' u (1000/$1):2 w l lt 3 lw 3 title "45 atm"

