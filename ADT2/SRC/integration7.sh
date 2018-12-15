#!/bin/bash

##############################################################################################################################################
#This script is written by Soumya Mukherjee, Bijit Mukherjee, Saikat Mukherjee, Subhankar Sardar and Satrajit Adhikari (corresponding author)
##############################################################################################################################################

M='NO OF STATE'

N1=`echo "$M*($M-1)/2" | bc -lq`

N=${N1%.*}

NS=`echo "$N+1" | bc -lq`

A1=`echo "10+3*$NS" | bc -lq`

A=${A1%.*}

B1=`echo "10+4*$NS" | bc -lq`

B=${B1%.*}

C1=`echo "10+5*$NS" | bc -lq`

C=${C1%.*}

D1=`echo "10+6*$NS" | bc -lq`

D=${D1%.*}

E1=`echo "2+$M" | bc -lq`

E=${E1%.*}

#parameters of 8th order Runge-Kutta methods
echo "      S=DSQRT(21.0D0)                                                   " >> integration.h
echo "*                                                                       " >> integration.h
echo "      CX1=0.5D0                                                         " >> integration.h 
echo "      CX2=(7.0D0+S)/14.0D0                                              " >> integration.h
echo "      CX3=(7.0D0-S)/14.0D0                                              " >> integration.h     
echo "*                                                                       " >> integration.h
echo "      CY21=0.5D0                                                        " >> integration.h 
echo "      CY31=0.25D0                                                       " >> integration.h  
echo "      CY32=0.25D0                                                       " >> integration.h
echo "      CY41=1.0D0/7.0D0                                                  " >> integration.h
echo "      CY42=(-7.0D0-3.0D0*S)/98.0D0                                      " >> integration.h
echo "      CY43=(21.0D0+5.0D0*S)/49.0D0                                      " >> integration.h  
echo "      CY51=(11.0D0+S)/84.0D0                                            " >> integration.h
echo "      CY52=(18.0D0+4.0D0*S)/63.0D0                                      " >> integration.h
echo "      CY53=(21.0D0-S)/252.0D0                                           " >> integration.h
echo "      CY61=(5.0D0+S)/48.0D0                                             " >> integration.h    
echo "      CY62=(9.0D0+S)/36.0D0                                             " >> integration.h
echo "      CY63=(231.0D0+14.0D0*S)/360.0D0                                   " >> integration.h
echo "      CY64=(63.0D0-7.0D0*S)/80.0D0                                      " >> integration.h
echo "      CY71=(10.0D0-S)/42.0D0                                            " >> integration.h
echo "      CY72=(-432.0D0+92.0D0*S)/315.0D0                                  " >> integration.h
echo "      CY73=(633.0D0-145.0D0*S)/90.0D0                                   " >> integration.h
echo "      CY74=(-504.0D0+115.0D0*S)/70.0D0                                  " >> integration.h
echo "      CY75=(63.0D0-13.0D0*S)/35.0D0                                     " >> integration.h
echo "      CY81=1.0D0/14.0D0                                                 " >> integration.h
echo "      CY82=(14.0D0-3.0D0*S)/126.0D0                                     " >> integration.h 
echo "      CY83=(13.0D0-3.0D0*S)/63.0D0                                      " >> integration.h
echo "      CY84=1.0D0/9.0D0                                                  " >> integration.h  
echo "      CY91=1.0D0/32.0D0                                                 " >> integration.h
echo "      CY92=(91.0D0-21.0D0*S)/576.0D0                                    " >> integration.h
echo "      CY93=11.0D0/72.0D0                                                " >> integration.h
echo "      CY94=(-385.0D0-75.0D0*S)/1152.0D0                                 " >> integration.h
echo "      CY95=(63.0D0+13.0D0*S)/128.0D0                                    " >> integration.h
echo "      CY101=1.0D0/14.0D0                                                " >> integration.h
echo "      CY102=1.0D0/9.0D0                                                 " >> integration.h
echo "      CY103=(-733.0D0-147.0D0*S)/2205.0D0                               " >> integration.h
echo "      CY104=(515.0D0+111.0D0*S)/504.0D0                                 " >> integration.h
echo "      CY105=(-51.0D0-11.0D0*S)/56.0D0                                   " >> integration.h
echo "      CY106=(132.0D0+28.0D0*S)/245.0D0                                  " >> integration.h
echo "      CY111=(-42.0D0+7.0D0*S)/18.0D0                                    " >> integration.h 
echo "      CY112=(-18.0D0+28.0D0*S)/45.0D0                                   " >> integration.h 
echo "      CY113=(-273.0D0-53.0D0*S)/72.0D0                                  " >> integration.h
echo "      CY114=(301.0D0+53.0D0*S)/72.0D0                                   " >> integration.h
echo "      CY115=(28.0D0-28.0D0*S)/45.0D0                                    " >> integration.h
echo "      CY116=(49.0D0-7.0D0*S)/18.0D0                                     " >> integration.h
echo "*                                                                       " >> integration.h

echo "      H=(RR(2)-RR(1))                                                   " >> integration.h

#solving of ADT equations along each positive step of first coordinate at highest value of second coordinate 
echo "      TR=RR1(1)+0.5D0*H                                                 " >> integration.h
echo "      TP=PH1(MM2+2)                                                     " >> integration.h                                                      
echo "      TOUTP=PH(MM2)                                                     " >> integration.h  
echo "*                                                                       " >> integration.h
echo "      DO I=1,NN                                                         " >> integration.h
echo "      Y(I)=0.0D0     ! Initial values of ADT angles.                    " >> integration.h
echo "      ENDDO                                                             " >> integration.h    
echo "*                                                                       " >> integration.h
echo "      DO I=1,MM1                                                        " >> integration.h
echo "         DO K=1,NN                                                      " >> integration.h                  
echo "            Y1(K)=Y(K)                                                  " >> integration.h
echo "         ENDDO                                                          " >> integration.h
echo "         CALL FEXR(TR,Y1,FR)                                            " >> integration.h
echo "         DO K=1,NN                                                      " >> integration.h
echo "            AK1(K)=FR(K)*H                                              " >> integration.h
echo "            Y2(K)=Y(K)+AK1(K)*CY21                                      " >> integration.h  
echo "         ENDDO                                                          " >> integration.h
echo "         CALL FEXR(TR+CX1*H,Y2,FR)                                      " >> integration.h 
echo "         DO K=1,NN                                                      " >> integration.h 
echo "            AK2(K)=FR(K)*H                                              " >> integration.h 
echo "            Y3(K)=Y(K)+AK1(K)*CY31+AK2(K)*CY32                          " >> integration.h
echo "         ENDDO                                                          " >> integration.h  
echo "         CALL FEXR(TR+CX1*H,Y3,FR)                                      " >> integration.h
echo "         DO K=1,NN                                                      " >> integration.h
echo "            AK3(K)=FR(K)*H                                              " >> integration.h
echo "            Y4(K)=Y(K)+AK1(K)*CY41+AK2(K)*CY42+AK3(K)*CY43              " >> integration.h 
echo "         ENDDO                                                          " >> integration.h
echo "         CALL FEXR(TR+CX2*H,Y4,FR)                                      " >> integration.h
echo "         DO K=1,NN                                                      " >> integration.h
echo "            AK4(K)=FR(K)*H                                              " >> integration.h    
echo "            Y5(K)=Y(K)+AK1(K)*CY51+AK3(K)*CY52+AK4(K)*CY53              " >> integration.h
echo "         ENDDO                                                          " >> integration.h
echo "         CALL FEXR(TR+CX2*H,Y5,FR)                                      " >> integration.h
echo "         DO K=1,NN                                                      " >> integration.h 
echo "            AK5(K)=FR(K)*H                                              " >> integration.h   
echo "            Y6(K)=Y(K)+AK1(K)*CY61+AK3(K)*CY62+AK4(K)*CY63+AK5(K)*CY64  " >> integration.h
echo "         ENDDO                                                          " >> integration.h  
echo "         CALL FEXR(TR+CX1*H,Y6,FR)                                      " >> integration.h
echo "         DO K=1,NN                                                      " >> integration.h  
echo "            AK6(K)=FR(K)*H                                              " >> integration.h
echo "            Y7(K)=Y(K)+AK1(K)*CY71+AK3(K)*CY72+AK4(K)*CY73+AK5(K)*CY74  " >> integration.h
echo "     $            +AK6(K)*CY75                                          " >> integration.h
echo "         ENDDO                                                          " >> integration.h
echo "         CALL FEXR(TR+CX3*H,Y7,FR)                                      " >> integration.h
echo "         DO K=1,NN                                                      " >> integration.h
echo "            AK7(K)=FR(K)*H                                              " >> integration.h
echo "            Y8(K)=Y(K)+AK1(K)*CY81+AK5(K)*CY82+AK6(K)*CY83+AK7(K)*CY84  " >> integration.h
echo "         ENDDO                                                          " >> integration.h
echo "         CALL FEXR(TR+CX3*H,Y8,FR)                                      " >> integration.h
echo "         DO K=1,NN                                                      " >> integration.h
echo "            AK8(K)=FR(K)*H                                              " >> integration.h 
echo "            Y9(K)=Y(K)+AK1(K)*CY91+AK5(K)*CY92+AK6(K)*CY93+AK7(K)*CY94  " >> integration.h
echo "     $            +AK8(K)*CY95                                          " >> integration.h      
echo "         ENDDO                                                          " >> integration.h
echo "         CALL FEXR(TR+CX1*H,Y9,FR)                                      " >> integration.h
echo "         DO K=1,NN                                                      " >> integration.h
echo "            AK9(K)=FR(K)*H                                              " >> integration.h
echo "            Y10(K)=Y(K)+AK1(K)*CY101+AK5(K)*CY102+AK6(K)*CY103          " >> integration.h
echo "     $             +AK7(K)*CY104+AK8(K)*CY105+AK9(K)*CY106              " >> integration.h  
echo "         ENDDO                                                          " >> integration.h
echo "         CALL FEXR(TR+CX2*H,Y10,FR)                                     " >> integration.h
echo "         DO K=1,NN                                                      " >> integration.h
echo "            AK10(K)=FR(K)*H                                             " >> integration.h
echo "            Y11(K)=Y(K)+AK5(K)*CY111+AK6(K)*CY112+AK7(K)*CY113          " >> integration.h
echo "     $             +AK8(K)*CY114+AK9(K)*CY115+AK10(K)*CY116             " >> integration.h 
echo "         ENDDO                                                          " >> integration.h
echo "         CALL FEXR(TR+H,Y11,FR)                                         " >> integration.h
echo "         DO K=1,NN                                                      " >> integration.h 
echo "            AK11(K)=FR(K)*H                                             " >> integration.h 
echo "         ENDDO                                                          " >> integration.h
echo "         DO K=1,NN                                                      " >> integration.h
echo "            Y(K)=Y(K)+(9.0D0*AK1(K)+49.0D0*AK8(K)+64.0D0*AK9(K)         " >> integration.h
echo "     $           +49.0D0*AK10(K)+9.0D0*AK11(K))/180.0D0                 " >> integration.h
echo "         ENDDO                                                          " >> integration.h
echo "         TR=RR(I)                                                       " >> integration.h  
echo "C                                                                       " >> integration.h
echo "         DO K=1,NN                                                      " >> integration.h
echo "           YY(I,K)=Y(K)                                                 " >> integration.h
echo "         ENDDO                                                          " >> integration.h
echo "C                                                                       " >> integration.h   
echo "      ENDDO                                                             " >> integration.h 
echo "CC                                                                      " >> integration.h
echo "CC                                                                      " >> integration.h

#solving of ADT equations along each negative step of second coordinate for each value of first coordinate
echo "      DO I=1,MM1                                                        " >> integration.h
echo "         TOUTR=RR(I)                                                    " >> integration.h
echo "         TR=TOUTR                                                       " >> integration.h
echo "         DO K=1,NN                                                      " >> integration.h
echo "            Y(K)=YY(I,K)                                                " >> integration.h
echo "         ENDDO                                                          " >> integration.h
echo "CCC                                                                     " >> integration.h
echo "         H=PH(1)-PH(2)                                                  " >> integration.h
echo "         TP=PH1(MM2+2)+0.5D0*H                                          " >> integration.h
echo "         DO J=MM2,1,-1                                                  " >> integration.h 
echo "           DO K=1,NN                                                    " >> integration.h
echo "              Y1(K)=Y(K)                                                " >> integration.h
echo "           ENDDO                                                        " >> integration.h
echo "           CALL FEXP(TP,Y1,FP)                                          " >> integration.h
echo "           DO K=1,NN                                                    " >> integration.h
echo "              AK1(K)=FP(K)*H                                            " >> integration.h
echo "              Y2(K)=Y(K)+AK1(K)*CY21                                    " >> integration.h
echo "           ENDDO                                                        " >> integration.h
echo "           CALL FEXP(TP+CX1*H,Y2,FP)                                    " >> integration.h
echo "           DO K=1,NN                                                    " >> integration.h
echo "              AK2(K)=FP(K)*H                                            " >> integration.h
echo "              Y3(K)=Y(K)+AK1(K)*CY31+AK2(K)*CY32                        " >> integration.h   
echo "           ENDDO                                                        " >> integration.h
echo "           CALL FEXP(TP+CX1*H,Y3,FP)                                    " >> integration.h
echo "           DO K=1,NN                                                    " >> integration.h
echo "              AK3(K)=FP(K)*H                                            " >> integration.h
echo "              Y4(K)=Y(K)+AK1(K)*CY41+AK2(K)*CY42+AK3(K)*CY43            " >> integration.h 
echo "           ENDDO                                                        " >> integration.h
echo "           CALL FEXP(TP+CX2*H,Y4,FP)                                    " >> integration.h
echo "           DO K=1,NN                                                    " >> integration.h
echo "              AK4(K)=FP(K)*H                                            " >> integration.h 
echo "              Y5(K)=Y(K)+AK1(K)*CY51+AK3(K)*CY52+AK4(K)*CY53            " >> integration.h
echo "           ENDDO                                                        " >> integration.h  
echo "           CALL FEXP(TP+CX2*H,Y5,FP)                                    " >> integration.h 
echo "           DO K=1,NN                                                    " >> integration.h
echo "              AK5(K)=FP(K)*H                                            " >> integration.h
echo "              Y6(K)=Y(K)+AK1(K)*CY61+AK3(K)*CY62+AK4(K)*CY63+AK5(K)*CY64" >> integration.h 
echo "           ENDDO                                                        " >> integration.h
echo "           CALL FEXP(TP+CX1*H,Y6,FP)                                    " >> integration.h 
echo "           DO K=1,NN                                                    " >> integration.h
echo "              AK6(K)=FP(K)*H                                            " >> integration.h
echo "              Y7(K)=Y(K)+AK1(K)*CY71+AK3(K)*CY72+AK4(K)*CY73+AK5(K)*CY74" >> integration.h 
echo "     $              +AK6(K)*CY75                                        " >> integration.h
echo "           ENDDO                                                        " >> integration.h
echo "           CALL FEXP(TP+CX3*H,Y7,FP)                                    " >> integration.h
echo "           DO K=1,NN                                                    " >> integration.h
echo "              AK7(K)=FP(K)*H                                            " >> integration.h
echo "              Y8(K)=Y(K)+AK1(K)*CY81+AK5(K)*CY82+AK6(K)*CY83+AK7(K)*CY84" >> integration.h 
echo "           ENDDO                                                        " >> integration.h 
echo "           CALL FEXP(TP+CX3*H,Y8,FP)                                    " >> integration.h 
echo "           DO K=1,NN                                                    " >> integration.h  
echo "              AK8(K)=FP(K)*H                                            " >> integration.h
echo "              Y9(K)=Y(K)+AK1(K)*CY91+AK5(K)*CY92+AK6(K)*CY93+AK7(K)*CY94" >> integration.h 
echo "     $              +AK8(K)*CY95                                        " >> integration.h
echo "           ENDDO                                                        " >> integration.h
echo "           CALL FEXP(TP+CX1*H,Y9,FP)                                    " >> integration.h
echo "           DO K=1,NN                                                    " >> integration.h
echo "              AK9(K)=FP(K)*H                                            " >> integration.h 
echo "              Y10(K)=Y(K)+AK1(K)*CY101+AK5(K)*CY102+AK6(K)*CY103        " >> integration.h
echo "     $               +AK7(K)*CY104+AK8(K)*CY105+AK9(K)*CY106            " >> integration.h  
echo "           ENDDO                                                        " >> integration.h 
echo "           CALL FEXP(TP+CX2*H,Y10,FP)                                   " >> integration.h
echo "           DO K=1,NN                                                    " >> integration.h
echo "              AK10(K)=FP(K)*H                                           " >> integration.h
echo "              Y11(K)=Y(K)+AK5(K)*CY111+AK6(K)*CY112+AK7(K)*CY113        " >> integration.h 
echo "     $               +AK8(K)*CY114+AK9(K)*CY115+AK10(K)*CY116           " >> integration.h
echo "           ENDDO                                                        " >> integration.h
echo "           CALL FEXP(TP+H,Y11,FP)                                       " >> integration.h
echo "           DO K=1,NN                                                    " >> integration.h
echo "              AK11(K)=FP(K)*H                                           " >> integration.h
echo "           ENDDO                                                        " >> integration.h
echo "           DO K=1,NN                                                    " >> integration.h
echo "              Y(K)=Y(K)+(9.0D0*AK1(K)+49.0D0*AK8(K)+64.0D0*AK9(K)       " >> integration.h  
echo "     $             +49.0D0*AK10(K)+9.0D0*AK11(K))/180.0D0               " >> integration.h
echo "           ENDDO                                                        " >> integration.h 
echo "           TP=PH(J)                                                     " >> integration.h 
echo "C                                                                       " >> integration.h
echo "           DO K=1,NN                                                    " >> integration.h
echo "             YYY(I,J,K)=Y(K)                                            " >> integration.h
echo "           ENDDO                                                        " >> integration.h
echo "*                                                                       " >> integration.h
echo "*                                                                       " >> integration.h

#printing of ADT angles in output files
echo "           DO K=1,NN                                                    " >> integration.h
echo "             WRITE($A+K,33)RR(I),PH(J),Y(K)                             " >> integration.h
echo "           ENDDO                                                        " >> integration.h
echo "C                                                                       " >> integration.h

#printing of ADT matrix elements in output files
echo "           CALL AMAT(Y,AA)                                              " >> integration.h
echo "           DO II=1,N                                                    " >> integration.h
echo "              WRITE($B+II,44)RR(I),PH(J),(AA(II,JJ),JJ=1,N)             " >> integration.h
echo "           ENDDO                                                        " >> integration.h
echo "*                                                                       " >> integration.h
echo "           DO II=1,N                                                    " >> integration.h
echo "              DO JJ=1,N                                                 " >> integration.h
echo "                IF(II.EQ.JJ)THEN                                        " >> integration.h
echo "                   UMAT(II,JJ)=UU(I,J,II)                               " >> integration.h
echo "                ELSE                                                    " >> integration.h
echo "                   UMAT(II,JJ)=0.0D0                                    " >> integration.h
echo "                ENDIF                                                   " >> integration.h
echo "              ENDDO                                                     " >> integration.h
echo "           ENDDO                                                        " >> integration.h
echo "*                                                                       " >> integration.h

#printing of diabatic potential energy matrix elements in output files
echo "           CALL WMAT(N,UMAT,AA,WA)                                      " >> integration.h
echo "           DO II=1,N                                                    " >> integration.h
echo "              WRITE($C+II,44)RR(I),PH(J),(WA(II,JJ),JJ=1,N)             " >> integration.h
echo "           ENDDO                                                        " >> integration.h
echo "*                                                                       " >> integration.h
echo "         ENDDO                                                          " >> integration.h
echo "                                                                        " >> integration.h
echo "         DO K=1,NN                                                      " >> integration.h
echo "           WRITE($A+K,*)                                                " >> integration.h
echo "         ENDDO                                                          " >> integration.h
echo "                                                                        " >> integration.h
echo "         DO II=1,N                                                      " >> integration.h
echo "           WRITE($B+II,*)                                               " >> integration.h
echo "           WRITE($C+II,*)                                               " >> integration.h
echo "         ENDDO                                                          " >> integration.h
echo "                                                                        " >> integration.h
echo "      ENDDO                                                             " >> integration.h
echo "CCC                                                                     " >> integration.h

#printing of residue of ADT angles in output files
echo "      DO K=1,NN                                                         " >> integration.h
echo "         DO I=1,MM1                                                     " >> integration.h                                                  
echo "            SS(K)=0.0D0                                                 " >> integration.h
echo "            DO J=1,MM2                                                  " >> integration.h
echo "               SS(K)=SS(K)+YYY(I,J,K)*DP                                " >> integration.h
echo "            ENDDO                                                       " >> integration.h
echo "            WRITE($D+K,55)RR(I),SS(K)                                   " >> integration.h   
echo "         ENDDO                                                          " >> integration.h
echo "      ENDDO                                                             " >> integration.h
echo "*                                                                       " >> integration.h
echo "   33 FORMAT(1X,3F12.6)                                                 " >> integration.h
echo "   44 FORMAT(1X,${E}F12.6)                                              " >> integration.h
echo "   55 FORMAT(1X,2F12.6)                                                 " >> integration.h
