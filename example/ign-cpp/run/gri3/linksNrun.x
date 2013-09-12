#!/bin/bash

RUNDIR=`pwd`
SRCDIR="$RUNDIR/../../src"
DATDIR="$RUNDIR/../../../../data"
PYDIR="$RUNDIR/../../../python"

echo "--------------------------------------------------------"
echo "run directory   : $RUNDIR"
echo "src directory   : $SRCDIR"
echo "python directory: $PYDIR"
echo "--------------------------------------------------------"
echo

# cd $SRCDIR; make clean; make output=yes; cd $RUNDIR
ln -fs $DATDIR/periodictable.dat .

if [ -f $DATDIR/chem_gri3.inp ]
then
ln -fs $DATDIR/chem_gri3.inp chem.inp
else
echo "Download grimech30.dat from:" 
echo "   http://www.me.berkeley.edu/gri_mech/version30/files30/grimech30.dat"
echo "then rename it to:"
echo "   $DATDIR/chem_gri3.inp"
exit
fi
if [ -f $DATDIR/therm_gri3.dat ]
then
ln -fs $DATDIR/therm_gri3.dat therm.dat
else
echo "Download thermo30.dat from:" 
echo "   http://www.me.berkeley.edu/gri_mech/version30/files30/thermo30.dat"
echo "then rename it to:"
echo "   $DATDIR/therm_gri3.inp"
exit
fi

ln -fs $SRCDIR/ign .

/bin/cp -f input.setup input.dat
./ign < input.dat

gnuplot plot.gp
./cleanfig.x specMaj.eps
./cleanfig.x specMin.eps
