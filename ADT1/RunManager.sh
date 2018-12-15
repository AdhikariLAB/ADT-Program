#!/bin/bash

##############################################################################################################################################
#This script is written by Soumya Mukherjee, Bijit Mukherjee, Saikat Mukherjee, Subhankar Sardar and Satrajit Adhikari (corresponding author)
##############################################################################################################################################

#################################
#Taking input from "Interface.sh"
#################################

NSTATE=$1
METHOD=$2

#####################################################
#Deleting previous directories and creating new ones
#####################################################

SrcDir=SRC
WrkDir=RUN
ResDir=RESULT

#deletion of 'RUN' directory from 'ADT1' folder, if exists
if [ -d $WrkDir ]; then

   rm -rf $WrkDir

fi

#deletion of 'RESULT' directory from 'ADT' folder, if exists
if [ -d ../$ResDir ]; then

   rm -rf ../$ResDir

fi

#creating new 'RUN' and 'RESULT' directories
mkdir -p $WrkDir
mkdir -p ../$ResDir

####################################################################################################
#Copying necessary files from source directory to working directory and performing necessary changes
####################################################################################################

#copying necessary python files from 'SRC' to 'RUN'
cp $SrcDir/adt$METHOD.py $WrkDir/
cp $SrcDir/ModuleBase.py $WrkDir/

#modification of the main python file according to user-defined number of states
sed -e "s/N = 'NO OF STATE'/N = $NSTATE/g" $WrkDir/adt$METHOD.py > aux 
mv aux $WrkDir/adt$METHOD.py

##########################
#Executing the python file
##########################

#entry into 'RUN' 
cd $WrkDir

#making the main python file executable and then, execution of that file
# chmod +x adt$METHOD.py
# ./adt$METHOD.py
python adt$METHOD.py $NSTATE

#exit from 'RUN'
cd ../

##################################################
#Transferring the output files to result directory
##################################################

#tranferring the output files to 'RESULT'
mv $WrkDir/*.DAT ../$ResDir
 

