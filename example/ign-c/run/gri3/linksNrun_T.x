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

#
#-Link executable
#
cd $SRCDIR; make clean; make output=yes gettig=yes; cd $RUNDIR
ln -fs $SRCDIR/ign .

#
#-Link data
#
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


# Set input file parameters
pfac=10
NiterMax=2000000
oFreq=1
deltat=1.e-3
deltatMax=1.0
tEnd=50.0 
Tini=1000.0
Temp_id=1500.0
deltaTemp=100.0
mech=chem.inp
thermo=therm.dat
withTab=0
getIgnDel=1
eqrat=1.0

# Lists of temperatures and pressures
Tinilist="800.0 810.0 820.0 830.0 840.0 850.0 860.0 870.0 880.0 890.0 900.0 910.0 920.0 930.0 940.0 950.0 975.0 1000.0 1025.0 1050.0 1075.0 1100.0 1125.0 1150.0 1175.0 1200.0 1250.0 1300.0"
pfaclist="2 10 50"

for ip in $pfaclist
do

touch igndel.dat; /bin/rm -f igndel.dat; touch igndel.dat
for it in $Tinilist
do
   echo "pfac          $ip"        >  input.dat
   echo "NiterMax      $NiterMax"  >> input.dat
   echo "oFreq         $oFreq"     >> input.dat
   echo "deltat        $deltat"    >> input.dat
   echo "deltatMax     $deltatMax" >> input.dat
   echo "tEnd          $tEnd"      >> input.dat
   echo "Tini          $it"        >> input.dat
   echo "Temp_id       $Temp_id"   >> input.dat
   echo "deltaTemp     $deltaTemp" >> input.dat
   echo "mech          $mech"      >> input.dat
   echo "thermo        $thermo"    >> input.dat
   echo "withTab       $withTab"   >> input.dat
   echo "getIgnDel     $getIgnDel" >> input.dat
   $PYDIR/pmix.py 1 $eqrat 1       >> input.dat
   echo "END"                      >> input.dat
   ./ign 
   echo "$it "> tinitmp.dat
   paste tinitmp.dat tid.dat >> igndel.dat
done

ignfile=igndel_${ip}_t.dat
if [ -e $ignfile ]
then
  mv $ignfile ${ignfile}.old
fi
mv igndel.dat $ignfile

done

gnuplot plot_id.gp
./cleanfig.x gri3_igndel.eps

exit

