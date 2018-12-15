#!/bin/bash

##############################################################################################################################################
#This script is written by Soumya Mukherjee, Bijit Mukherjee, Saikat Mukherjee, Subhankar Sardar and Satrajit Adhikari (corresponding author)
##############################################################################################################################################

#################################
#Taking input from "Interface.sh"
#################################

NSTATE=$1
RLIM=$2
PLIM=$3
INTEGRATION=$4
COMPILER=$5

##############################################################
#Deleting previous files and directories and creating new ones
##############################################################

InpDir=INPUT
SrcDir=SRC
WrkDir=RUN
ResDir=RESULT

#deletion of compile script, 'MakeCompile.sh', if exists
if [ -f MakeCompile.sh ]; then

   rm MakeCompile.sh

fi

#deletion of 'RUN' directory from 'ADT2' folder, if exists
if [ -d $WrkDir ]; then

   rm -rf $WrkDir

fi

#deletion of 'RESULT' directory from 'ADT' folder, if exists
if [ -d ../$ResDir ]; then

   rm -rf ../$ResDir

fi

#deletion of any existing executable file created from previous runs
for f in *.out; do

   [ -e "$f" ] && rm $f

done

###############################################################
#Calculating the number of possible nonadiabatic coupling terms
###############################################################

N1=`echo "$NSTATE*($NSTATE-1)/2" | bc -lq`

N=${N1%.*}

########################################################################################################################
#Copying necessary files from source directory and input directory to working directory and performing necessary changes
########################################################################################################################

#creating new 'RUN' and 'RESULT' directories
mkdir -p $WrkDir
mkdir -p ../$ResDir

#making the shell script files executable 
chmod +x $SrcDir/*.sh

#copying necessary fortran and bash script files from 'SRC' to 'RUN'
cp $SrcDir/adt.f $WrkDir/
cp $SrcDir/open.sh $WrkDir/
cp $SrcDir/integration$INTEGRATION.sh $WrkDir/

#copying necessary input files from 'INPUT' to 'RUN'
for i in `seq 1 1 $NSTATE`
    do

      cp ../$InpDir/ADIABATICPES_$i.DAT $WrkDir/ 

    done

for i in `seq 1 1 $N`
    do

      cp ../$InpDir/TAUR_$i.DAT $WrkDir/ 
      cp ../$InpDir/TAUP_$i.DAT $WrkDir/

    done

#modifications of the main fortran file according to user-defined number of states, 
#number of grid points and integration path
sed -e "s/NN='NO OF TAU'/NN=$N/g" $WrkDir/adt.f > aux
sed -e "s/N='NO OF STATE'/N=$NSTATE/g" aux > $WrkDir/adt.f

sed -e "s/NNN='NO OF TAU'/NNN=$N/g" $WrkDir/adt.f > aux
mv aux $WrkDir/adt.f

sed -e "s/MM1='C1 LIMIT'/MM1=$RLIM/g" $WrkDir/adt.f > aux
sed -e "s/MM2='C2 LIMIT'/MM2=$PLIM/g" aux > $WrkDir/adt.f

sed -e "s/NSTATE = '\$NSTATE\$'/NSTATE = $NSTATE/g" $WrkDir/adt.f > aux 
mv aux $WrkDir/adt.f

if [ $INTEGRATION -le 4 ]
   
   then 

   sed -e "s/YY('NO OF INITIAL GRID',NN)/YY(MM2,NN)/g" $WrkDir/adt.f > aux 
   mv aux $WrkDir/adt.f

else

   sed -e "s/YY('NO OF INITIAL GRID',NN)/YY(MM1,NN)/g" $WrkDir/adt.f > aux 
   mv aux $WrkDir/adt.f

fi

#modification of 'open.sh' according to user-defined number of states, 
sed -e "s/M='NO OF STATE'/M=$NSTATE/g" $WrkDir/open.sh > aux
mv aux $WrkDir/open.sh

#modification of the relevant 'integration' script according to user-defined 
#number of states 
sed -e "s/M='NO OF STATE'/M=$NSTATE/g" $WrkDir/integration$INTEGRATION.sh > aux
mv aux $WrkDir/integration$INTEGRATION.sh

##########################################################################
#Running the script files and fortran file to generate the executable file
##########################################################################

#entry into 'RUN' 
cd $WrkDir

#making 'open.sh' executable and then, execution of that file
chmod +x open.sh
./open.sh

#making relevant 'integration' script file executable and then, execution of that file
chmod +x integration$INTEGRATION.sh
./integration$INTEGRATION.sh

#exit from 'RUN'
cd ../

#creating a compilation script, 'MakeCompile.sh'
echo -e "#!/bin/bash\n" >> MakeCompile.sh

if [ $COMPILER -eq 1 ]

   then

   echo -n "ifort $WrkDir/adt.f -o Execute"$NSTATE"st.out" >> MakeCompile.sh

elif [ $COMPILER -eq 2 ]

   then

   echo -n "gfortran $WrkDir/adt.f -o Execute"$NSTATE"st.out" >> MakeCompile.sh

fi

#making 'MakeCompile.sh' executable and then, execution of that file
chmod +x MakeCompile.sh
./MakeCompile.sh

#deletion of 'MakeCompile.sh'
rm MakeCompile.sh 

#################################
#Running ./Execute"$NSTATE"st.out
#################################

mv Execute"$NSTATE"st.out $WrkDir

cd $WrkDir

./Execute"$NSTATE"st.out

cd ../

##################################################
#Transferring the output files to result directory
##################################################

#tranferring the output files to 'RESULT'
for i in `seq 1 1 $N`
   do 
  
     mv $WrkDir/ADT_ANGLE_$i.DAT ../$ResDir
 
   done

for i in `seq 1 1 $NSTATE`
   do 
 
     mv $WrkDir/AMAT_ROW_$i.DAT ../$ResDir/
     mv $WrkDir/DIAPES_$i.DAT ../$ResDir/

   done

for i in `seq 1 1 $N`
   do 
 
     mv $WrkDir/ADT_ANGLE_RES_$i.DAT ../$ResDir

   done

     mv $WrkDir/INFORMATION.DAT ../$ResDir

