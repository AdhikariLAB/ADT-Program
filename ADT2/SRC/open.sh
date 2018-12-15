#!/bin/bash

##############################################################################################################################################
#This script is written by Soumya Mukherjee, Bijit Mukherjee, Saikat Mukherjee, Subhankar Sardar and Satrajit Adhikari (corresponding author)
##############################################################################################################################################

M='NO OF STATE'

#calculating number of NACTs
N1=`echo "$M*($M-1)/2" | bc -lq`

N=${N1%.*}

NS=`echo "$N+1" | bc -lq`

#deletion of header file, 'openfile.h', if exists
if [ -f openfile.h ]; then

rm openfile.h

fi

##################################################
#Printing of 'OPEN' statements in 'openfile.h'
##################################################

#opening of adiabatic energy files
for i in `seq 1 1 $M`
do

A1=`echo "10+$i" | bc -lq`

A=${A1%.*}

echo -e "       OPEN($A,FILE='ADIABATICPES_$i.DAT',STATUS='OLD')" >> openfile.h

done

echo  >> openfile.h

#opening of NACT (along first coordinate) files
for i in `seq 1 1 $N`
do 

B1=`echo "10+$NS+$i" | bc -lq`

B=${B1%.*}

echo -e "       OPEN($B,FILE='TAUR_$i.DAT',STATUS='OLD')" >> openfile.h

done

echo  >> openfile.h

#opening of NACT (along second coordinate) files
for i in `seq 1 1 $N`
do 

C1=`echo "10+2*$NS+$i" | bc -lq`

C=${C1%.*}

echo -e "       OPEN($C,FILE='TAUP_$i.DAT',STATUS='OLD')" >> openfile.h

done

echo  >> openfile.h

#opening of ADT angle files
for i in `seq 1 1 $N`
do

D1=`echo "10+3*$NS+$i" | bc -lq`

D=${D1%.*}

echo -e "       OPEN($D,FILE='ADT_ANGLE_$i.DAT',STATUS='UNKNOWN')" >> openfile.h

done

echo  >> openfile.h

#opening of ADT matrix files
for i in `seq 1 1 $M`
do

E1=`echo "10+4*$NS+$i" | bc -lq`

E=${E1%.*}

echo -e "       OPEN($E,FILE='AMAT_ROW_$i.DAT',STATUS='UNKNOWN')" >> openfile.h

done

echo  >> openfile.h

#opening of diabatic potential energy matrix files
for i in `seq 1 1 $M`
do

F1=`echo "10+5*$NS+$i" | bc -lq`

F=${F1%.*}

echo -e "       OPEN($F,FILE='DIAPES_$i.DAT',STATUS='UNKNOWN')" >> openfile.h

done

echo  >> openfile.h

#opening of ADT angle residue files
for i in `seq 1 1 $N`
do

G1=`echo "10+6*$NS+$i" | bc -lq`

G=${G1%.*}

echo -e "       OPEN($G,FILE='ADT_ANGLE_RES_$i.DAT',STATUS='UNKNOWN')" >> openfile.h

done

echo  >> openfile.h

############################################
#Reading of input data from input files
############################################

H1=`echo "10+$NS" | bc -lq`

H=${H1%.*}

I1=`echo "10+2*$NS" | bc -lq`

I=${I1%.*}

echo -e "       DO I=1,MM1"  >> openfile.h

echo -e "          DO J=1,MM2" >> openfile.h
 
echo -e "             DO K=1,N" >> openfile.h

echo -e "                READ(10+K,*)RR(I),PH(J),UU(I,J,K)" >> openfile.h

echo -e "             ENDDO" >> openfile.h

echo -e "             DO K=1,NN" >> openfile.h

echo -e "                READ($H+K,*)R,P,TAUR(I,J,K)" >> openfile.h

echo -e "                READ($I+K,*)R,P,TAUP(I,J,K)" >> openfile.h

echo -e "             ENDDO" >> openfile.h

echo -e "          ENDDO" >> openfile.h

echo -e "          DO K=1,N" >> openfile.h

echo -e "            READ(10+K,*)" >> openfile.h

echo -e "          ENDDO" >> openfile.h

echo -e "          DO K=1,NN" >> openfile.h

echo -e "            READ($H+K,*)" >> openfile.h

echo -e "            READ($I+K,*)" >> openfile.h

echo -e "          ENDDO" >> openfile.h

echo -e "       ENDDO" >> openfile.h

####################################################################
#Printing of authors' information and description of output data
####################################################################

echo -e                                                                                                              >> INFORMATION.DAT

echo -e "################################################################################################################################ " \
>> INFORMATION.DAT

echo -e "#Authors are Soumya Mukherjee, Bijit Mukherjee, Saikat Mukherjee, Subhankar Sardar and Satrajit Adhikari (corresponding author)" \
>> INFORMATION.DAT

echo -e "################################################################################################################################ " \
>> INFORMATION.DAT

echo -e                                                                                                              >> INFORMATION.DAT

echo -e "############################################################ "                                              >> INFORMATION.DAT

echo -e "#This file provides a brief description of the result files "                                               >> INFORMATION.DAT

echo -e "############################################################# "                                             >> INFORMATION.DAT

echo -e                                                                                                              >> INFORMATION.DAT

#printing of information regarding the ADT angle files
echo -e "################# "                                                                                         >> INFORMATION.DAT

echo -e "#ADT Angle Files "                                                                                          >> INFORMATION.DAT

echo -e "################# "                                                                                         >> INFORMATION.DAT

echo -e                                                                                                              >> INFORMATION.DAT

echo -e "----------------------------------------------------------------------------------------------------------" >> INFORMATION.DAT

echo -e "        File Name               First Column              Second Column           Third Column            " >> INFORMATION.DAT

echo -e "----------------------------------------------------------------------------------------------------------" >> INFORMATION.DAT

echo -e                                                                                                              >> INFORMATION.DAT

angle=0

for i in `seq 2 1 $M`
do

k=`echo "$i-1" |bc -lq`

k1=${k%.*}

for j in `seq 1 1 $k1`
do

angle=`echo "$angle+1" |bc -lq`

ang=${angle%.*}

echo -e "     ADT_ANGLE_$ang.DAT         first coordinate            second coordinate          theta($j,$i)   "     >> INFORMATION.DAT 

echo -e                                                                                                              >> INFORMATION.DAT

done

done

echo -e "----------------------------------------------------------------------------------------------------------" >> INFORMATION.DAT

echo -e                                                                                                              >> INFORMATION.DAT

#printing of information regarding the ADT matrix files
echo -e "################## "                                                                                        >> INFORMATION.DAT

echo -e "#ADT Matrix Files "                                                                                         >> INFORMATION.DAT

echo -e "################## "                                                                                        >> INFORMATION.DAT

echo -e                                                                                                              >> INFORMATION.DAT

echo -e "----------------------------------------------------------------------------------------------------------" >> INFORMATION.DAT

echo -e "        File Name               First Column              Second Column           Other Columns           " >> INFORMATION.DAT

echo -e "----------------------------------------------------------------------------------------------------------" >> INFORMATION.DAT

echo -e                                                                                                              >> INFORMATION.DAT

angle=0

for i in `seq 1 1 $M`
do

echo -e "       AMAT_ROW_$i.DAT         first coordinate           second coordinate       A($i,1) to A($i,$M)   "   >> INFORMATION.DAT 

echo -e                                                                                                              >> INFORMATION.DAT

done

echo -e "----------------------------------------------------------------------------------------------------------" >> INFORMATION.DAT

echo -e                                                                                                              >> INFORMATION.DAT

#printing of information regarding the diabatic potential energy matrix files
echo -e "####################### "                                                                                   >> INFORMATION.DAT

echo -e "#Diabatic Matrix Files "                                                                                    >> INFORMATION.DAT

echo -e "####################### "                                                                                   >> INFORMATION.DAT

echo -e                                                                                                              >> INFORMATION.DAT

echo -e "----------------------------------------------------------------------------------------------------------" >> INFORMATION.DAT

echo -e "        File Name               First Column              Second Column           Other Columns           " >> INFORMATION.DAT

echo -e "----------------------------------------------------------------------------------------------------------" >> INFORMATION.DAT

echo -e                                                                                                              >> INFORMATION.DAT

angle=0

for i in `seq 1 1 $M`
do

echo -e "       DIAPES_$i.DAT           first coordinate           second coordinate      W($i,1) to W($i,$M)    "   >> INFORMATION.DAT 

echo -e                                                                                                              >> INFORMATION.DAT

done

echo -e "----------------------------------------------------------------------------------------------------------" >> INFORMATION.DAT

echo -e                                                                                                              >> INFORMATION.DAT

#printing of information regarding the ADT angle residue files
echo -e "######################### "                                                                                 >> INFORMATION.DAT

echo -e "#ADT Angle Residue Files "                                                                                  >> INFORMATION.DAT

echo -e "######################### "                                                                                 >> INFORMATION.DAT

echo -e                                                                                                              >> INFORMATION.DAT

echo -e "----------------------------------------------------------------------------------------------------------" >> INFORMATION.DAT

echo -e "        File Name               First Column              Second Column                                   " >> INFORMATION.DAT

echo -e "----------------------------------------------------------------------------------------------------------" >> INFORMATION.DAT

echo -e                                                                                                              >> INFORMATION.DAT

angle=0

for i in `seq 2 1 $M`
do

k=`echo "$i-1" |bc -lq`

k1=${k%.*}

for j in `seq 1 1 $k1`
do

angle=`echo "$angle+1" |bc -lq`

ang=${angle%.*}

echo -e "     ADT_ANGLE_RES_$ang.DAT      first coordinate            residue($j,$i)                            "     >> INFORMATION.DAT 

echo -e                                                                                                              >> INFORMATION.DAT

done

done

echo -e "----------------------------------------------------------------------------------------------------------" >> INFORMATION.DAT


