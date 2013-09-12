#!/bin/bash

gnupl=/usr/local/gnuplot-4.4.2/bin/gnuplot
$gnupl ../plot_id1.gp
$gnupl ../plot_id2.gp
../cleanfig.x iOct_id1.eps
../cleanfig.x iOct_id2.eps

exit

