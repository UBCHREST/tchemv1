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
#cd $SRCDIR; make clean; make gettig=yes; cd $RUNDIR
ln -fs $SRCDIR/ign .

#
#-Link data
#
ln -fs $DATDIR/periodictable.dat .
if [ -f $DATDIR/chem_isoOct.inp ]
then
ln -fs $DATDIR/chem_isoOct.inp chem.inp
else
echo "Download ic8_ver3_mech.txt from:" 
echo "   https://www-pls.llnl.gov/data/docs/science_and_technology/chemistry/combustion/ic8_ver3_mech.txt"
echo "then rename it to:"
echo "   $DATDIR/chem_isoOct.inp"
exit
fi
if [ -f $DATDIR/therm_isoOct.dat ]
then
ln -fs $DATDIR/therm_isoOct.dat therm.dat
else
echo "Download prf_v3_therm_dat.txt from:" 
echo "   https://www-pls.llnl.gov/data/docs/science_and_technology/chemistry/combustion/prf_v3_therm_dat.txt"
echo "then rename it to:"
echo "   $DATDIR/therm_isoOct.inp"
exit
fi

#
#-Set input file parameters
#
pfac=10
NiterMax=2000000
oFreq=1
deltat=10.0
deltatMax=10.0
tEnd=10.0 
Tini=1000.0
Temp_id=1500.0
deltaTemp=10000.0
mech=chem.inp
thermo=therm.dat
withTab=0
getIgnDel=1

#
#-Lists of temperatures and pressures
#
Tinilist="650.0 660.0 670.0 680.0 690.0 700.0 710.0 720.0 730.0 740.0 750.0 760.0 770.0 780.0 790.0 800.0 810.0 820.0 830.0 840.0 850.0 860.0 870.0 880.0 890.0 900.0 910.0 920.0 930.0 940.0 950.0 975.0 1000.0 1025.0 1050.0 1075.0 1100.0 1125.0 1150.0 1175.0 1200.0 1250.0 1300.0"
pfaclist="13 16 34 40 45"

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
   echo "spec  IC8H18  0.0164837326409 " >> input.dat
   echo "spec  O2      0.206046658012  " >> input.dat
   echo "spec  N2      0.767929501554  " >> input.dat
   echo "spec  AR      0.00954010779338" >> input.dat
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
./cleanfig.x igndel_t.eps

exit

