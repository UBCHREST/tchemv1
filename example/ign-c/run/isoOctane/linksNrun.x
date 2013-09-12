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


cd $SRCDIR; make clean; make output=yes; cd $RUNDIR
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

ln -fs $SRCDIR/ign .

# Set input file parameters
pfac=45
NiterMax=2000000
oFreq=10
deltat=1.e-10
deltatMax=1.e-3
tEnd=2.0 
Tini=1000.0
Temp_id=1500.0
deltaTemp=1.0
mech=chem.inp
thermo=therm.dat
withTab=0
getIgnDel=1
# list of equivalence ratios
eqrat="0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.2 2.4 2.6 2.8 3.0"
#eqrat="0.3 0.4"

touch igndel.dat; /bin/rm -f igndel.dat; touch igndel.dat
for i in $eqrat
do
   echo "pfac          $pfac"      >  input.dat
   echo "NiterMax      $NiterMax"  >> input.dat
   echo "oFreq         $oFreq"     >> input.dat
   echo "deltat        $deltat"    >> input.dat
   echo "deltatMax     $deltatMax" >> input.dat
   echo "tEnd          $tEnd"      >> input.dat
   echo "Tini          $Tini"      >> input.dat
   echo "Temp_id       $Temp_id"   >> input.dat
   echo "deltaTemp     $deltaTemp" >> input.dat
   echo "mech          $mech"      >> input.dat
   echo "thermo        $thermo"    >> input.dat
   echo "withTab       $withTab"   >> input.dat
   echo "getIgnDel     $getIgnDel" >> input.dat
   $PYDIR/pmix.py 8 $i 1           >> input.dat
   echo "END"                      >> input.dat
   ./ign
   echo "$i " > eqtmp.dat
   paste eqtmp.dat tid.dat >> igndel.dat
done

ignfile=igndel_${pfac}atm.dat
if [ -e $ignfile ]
then
  mv $ignfile ${ignfile}.old
fi
mv igndel.dat $ignfile

exit
