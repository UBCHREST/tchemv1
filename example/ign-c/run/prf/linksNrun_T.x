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
cd $SRCDIR; make clean; make output=yes; cd $RUNDIR
ln -fs $SRCDIR/ign .

#
#-Link data
#
ln -fs $DATDIR/periodictable.dat .
if [ -f $DATDIR/chem_prf.inp ]
then
ln -fs $DATDIR/chem_prf.inp chem.inp
else
echo "Download prf_2d_mech.txt from:" 
echo "   https://www-pls.llnl.gov/data/docs/science_and_technology/chemistry/combustion/prf_2d_mech.txt"
echo "then rename it to:"
echo "   $DATDIR/chem_prf.inp"
exit
fi
if [ -f $DATDIR/therm_prf.dat ]
then
ln -fs $DATDIR/therm_prf.dat therm.dat
else
echo "Download prf_2d_therm.txt from:" 
echo "   https://www-pls.llnl.gov/data/docs/science_and_technology/chemistry/combustion/prf_2d_therm.txt"
echo "then rename it to:"
echo "   $DATDIR/therm_prf.inp"
exit
fi

#
#-Set input file parameters
#
pfac=1.0
NiterMax=2000000
oFreq=1
deltat=10.0
deltatMax=10.0
tEnd=10.0 
Tini=1000.0
Temp_id=1500.0
deltaTemp=1000.0
mech=chem.inp
thermo=therm.dat
withTab=0
getIgnDel=1

#
#-Lists of temperatures, pressures, and octane numbers
#
Tinilist="650.0 660.0 670.0 680.0 690.0 700.0 710.0 720.0 730.0 740.0 750.0 760.0 770.0 780.0 790.0 800.0 810.0 820.0 830.0 840.0 850.0 860.0 870.0 880.0 890.0 900.0 910.0 920.0 930.0 940.0 950.0 975.0 1000.0 1025.0 1050.0 1075.0 1100.0 1125.0 1150.0 1175.0 1200.0 1250.0 1300.0"
pfaclist="15 30 45"
onlist="100 90 80 70 60 50"

for ion in $onlist
do

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
   ${PYDIR}/pmixPRF.py $ion 1      >> input.dat
   echo "END"                      >> input.dat
   ./ign 
   echo "$it "> tinitmp.dat
   paste tinitmp.dat tid.dat >> igndel.dat
done # loop over temperature list

ignfile=igndel_p${ip}atm_on${ion}_t.dat
if [ -e $ignfile ]
then
  mv $ignfile ${ignfile}.old
fi
mv igndel.dat $ignfile

done # loop over pressure list

done # loop over the octane number list

exit

gnuplot plot_id.gp
./cleanfig.x igndel_t.eps

exit

