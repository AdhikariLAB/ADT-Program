#!/bin/bash

##############################################################################################################################################
#This script is written by Soumya Mukherjee, Bijit Mukherjee, Saikat Mukherjee, Subhankar Sardar and Satrajit Adhikari (corresponding author)
##############################################################################################################################################

initial=$(date +%s.%N)

###########################
#Taking input from the user
###########################

#INPUT=$1
JOB=$1
NSTATE=$2
METHOD=$3
RLIM=$4
PLIM=$5
INTEGRATION=$6
COMPILER=$7

###############################
#Executing the necessary script
###############################

#running molpro input
#if [ $INPUT -eq 1 ]; then
#   molpro $INPUT
#fi

#making the 'ADT1'/'ADT2' directory editable
chmod +w ADT$JOB

#entry into 'ADT1'/'ADT2'
cd ADT$JOB

#Executing the script-file, 'RunManager.sh' inside 'ADT1'/'ADT2' directory
if [ $JOB -eq 1 ]; then

   ./RunManager.sh $NSTATE $METHOD

else

   ./RunManager.sh $NSTATE $RLIM $PLIM $INTEGRATION $COMPILER

fi

#exit from 'ADT1'/'ADT2'
cd ../

#making the 'ADT1'/'ADT2' directory non-editable
chmod -w ADT$JOB

#calculation and printing of CPU time
final=$(date +%s.%N)
dt=$(echo "$final - $initial" | bc)
dd=$(echo "$dt/86400" | bc)
dt2=$(echo "$dt-86400*$dd" | bc)
dh=$(echo "$dt2/3600" | bc)
dt3=$(echo "$dt2-3600*$dh" | bc)
dm=$(echo "$dt3/60" | bc)
ds=$(echo "$dt3-60*$dm" | bc)

echo "Total runtime  $dd days $dh hours $dm minutes $ds seconds" > time.log

