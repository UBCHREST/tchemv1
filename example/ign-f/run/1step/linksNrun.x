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


cd $SRCDIR; make exec1step; cd $RUNDIR

ln -fs $DATDIR/periodictable.dat .
ln -fs $DATDIR/chem_1step.inp chem.inp
ln -fs $DATDIR/therm_1step.dat therm.dat

/bin/cp input.setup input.dat
ln -fs $SRCDIR/ign-mod .

./ign-mod < input.dat
