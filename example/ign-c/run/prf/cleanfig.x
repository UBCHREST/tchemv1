#!/bin/bash
if [ $# -ne 1 ]
then
  echo "usage: $0 filename"
  exit 
fi
FIGFILE=$1
awk '{if ($1=="/LT0") print "/LT0 { PL [] 1 0 0 DL } def"; else print $0;}' $FIGFILE > tmpfig; mv tmpfig $FIGFILE
awk '{if ($1=="/LT1") print "/LT1 { PL [] 0 1 0 DL } def"; else print $0;}' $FIGFILE > tmpfig; mv tmpfig $FIGFILE
awk '{if ($1=="/LT2") print "/LT2 { PL [] 0 0 1 DL } def"; else print $0;}' $FIGFILE > tmpfig; mv tmpfig $FIGFILE
awk '{if ($1=="/LT3") print "/LT3 { PL [] 0 1 1 DL } def"; else print $0;}' $FIGFILE > tmpfig; mv tmpfig $FIGFILE
awk '{if ($1=="/LT4") print "/LT4 { PL [] 0 1 1 DL } def"; else print $0;}' $FIGFILE > tmpfig; mv tmpfig $FIGFILE
awk '{if ($1=="/LT5") print "/LT5 { PL [] 1 1 0 DL } def"; else print $0;}' $FIGFILE > tmpfig; mv tmpfig $FIGFILE

