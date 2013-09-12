#!/bin/bash

gnupl=/usr/local/gnuplot-4.4.2/bin/gnuplot
$gnupl ../plot_id.gp
../cleanfig.x prf_id1.eps
../cleanfig.x prf_id2.eps

exit

