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
ln -fs $DATDIR/chem_gri3.inp chem.inp
ln -fs $DATDIR/therm_gri3.dat therm.dat

cd $SRCDIR; make execgri3 output=yes; cd $RUNDIR
ln -fs $SRCDIR/ign .

./ign < input.dat

gnuplot plot.gp
./cleanfig.x specMaj.eps
./cleanfig.x specMin.eps

