
*****************************************************************************************************************************************************
C This fortran code is written by Soumya Mukherjee, Bijit Mukherjee, Saikat Mukherjee, Subhankar Sardar and Satrajit Adhikari (corresponding author)
*****************************************************************************************************************************************************

*************************************************MAIN-PROGRAM**************************************

********************************
C     Declaration of variables  
******************************** 

      IMPLICIT NONE
      INTEGER MM1,MM2                                                
      PARAMETER (MM1='C1 LIMIT')      !No. of first coordinate grid
      PARAMETER (MM2='C2 LIMIT')      !No. of second coordinate grid
      INTEGER MM,NN,N,NNN             !Total grid points, no. of ADT angles, no. of states, and no. of ADT angles in 'FEXP' and 'FEXR'
      INTEGER I,J,K                   !Loop index
      INTEGER II,JJ                   !Loop index
      PARAMETER (MM=MM1*MM2)
      PARAMETER (NN='NO OF TAU')              
      PARAMETER (N='NO OF STATE')                 
      PARAMETER (NNN='NO OF TAU')
      REAL*8 RR(MM1),PH(MM2),R,P    !First coordinate grid, second coordinate grid, temporary storage of first and second coordiante 
      REAL*8 DR,DP,H                !Interval in first and second coordinate
      REAL*8 TAUR(MM1,MM2,NN)       !Components of NACTs along first coordiante
      REAL*8 TAUP(MM1,MM2,NN)       !Components of NACTs along second coordinate
      REAL*8 UU(MM1,MM2,N)          !Adiabatic potential energy data for each grid point
      REAL*8 UMAT(N,N)              !Adiabatic potential energy in matrix form for a particular grid point
      REAL*8 SS(NN),Y(NN),YY('NO OF INITIAL GRID',NN),YYY(MM1,MM2,NN)   !ADT angle residue, ADT angles 
      REAL*8 Y1(NN),Y2(NN),Y3(NN),Y4(NN),Y5(NN),Y6(NN)                  !parameters of 8th order Runge-Kutta method 
      REAL*8 Y7(NN),Y8(NN),Y9(NN),Y10(NN),Y11(NN)                       !parameters of 8th order Runge-Kutta method 
      REAL*8 AK1(NN),AK2(NN),AK3(NN),AK4(NN),AK5(NN),AK6(NN)            !parameters of 8th order Runge-Kutta method 
      REAL*8 AK7(NN),AK8(NN),AK9(NN),AK10(NN),AK11(NN)                  !parameters of 8th order Runge-Kutta method
      REAL*8 FR(NN),FP(NN)                                              !parameters of 8th order Runge-Kutta method 
      REAL*8 S,CX1,CX2,CX3,CY21,CY31,CY32,CY41,CY42,CY43                !parameters of 8th order Runge-Kutta method  
      REAL*8 CY51,CY52,CY53,CY61,CY62,CY63,CY64                         !parameters of 8th order Runge-Kutta method
      REAL*8 CY71,CY72,CY73,CY74,CY75,CY81,CY82,CY83,CY84               !parameters of 8th order Runge-Kutta method
      REAL*8 CY91,CY92,CY93,CY94,CY95                                   !parameters of 8th order Runge-Kutta method
      REAL*8 CY101,CY102,CY103,CY104,CY105,CY106                        !parameters of 8th order Runge-Kutta method
      REAL*8 CY111,CY112,CY113,CY114,CY115,CY116                        !parameters of 8th order Runge-Kutta method
      REAL*8 AA(N,N),WA(N,N)                        !A matrices, diabatic potential energy matrices
      REAL*8 A(NN,N,N)                        !A matrices, diabatic potential energy matrices
      REAL*8 RR1,PH1                    !First and second coordinate grid after grid expansion
      REAL*8 TR,TP,TOUTR,TOUTP          !Variables related to first and second coordinates
      REAL*8 TAUR_1                     !Components of NACTs along first coordinate (used in 'COMMON' block)
      REAL*8 TAUP_1                     !Components of NACTs along second coordiante (used in 'COMMON' block)
      COMMON/GRID/RR1(MM1+2),PH1(MM2+2) !Common block for the two expanded coordinates
      COMMON/TRR/TAUR_1(MM1+2,MM2+2,NNN)!Common block for components of NACTs along first coordinate (grid expanded)
      COMMON/TPP/TAUP_1(MM1+2,MM2+2,NNN)!Common block for components of NACTs along second coordinate (grid expanded)
      COMMON/RHH/TOUTR                  !Common block for TOUTR
      COMMON/PHH/TOUTP                  !Common block for TOUTP

********************************
C     Opening Statements
********************************

      include 'openfile.h'

*********************************************************
C     Calculation of intervals of the two coordiantes
*********************************************************

      DR=RR(2)-RR(1)
      DP=PH(2)-PH(1)

************************************************************
C     Expansion of the grid and NACT components
************************************************************
 
      DO I=2,MM1+1
         RR1(I)=RR(I-1)
      ENDDO
      RR1(1)=RR1(2)-DR
      RR1(MM1+2)=RR1(MM1+1)+DR
      DO I=2,MM2+1
         PH1(I)=PH(I-1)
      ENDDO
      PH1(1)=PH1(2)-DP
      PH1(MM2+2)=PH1(MM2+1)+DP

      DO K=1,NN
        DO I=2,MM1+1
          DO J=2,MM2+1
            TAUR_1(I,J,K)=TAUR(I-1,J-1,K)
            TAUP_1(I,J,K)=TAUP(I-1,J-1,K)
            ENDDO
         ENDDO
      ENDDO

      DO K=1,NN
         TAUR_1(1,2,K)=TAUR_1(2,2,K)
         TAUR_1(MM1+2,2,K)=TAUR_1(MM1+1,2,K)
         TAUR_1(2,MM2+2,K)=TAUR_1(2,MM2+1,K)
         TAUR_1(MM1+1,MM2+2,K)=TAUR_1(MM1+1,MM2+1,K)
         TAUP_1(1,2,K)=TAUP_1(2,2,K)
         TAUP_1(MM1+2,2,K)=TAUP_1(MM1+1,2,K)
         TAUP_1(2,MM2+2,K)=TAUP_1(2,MM2+1,K)
         TAUP_1(MM1+1,MM2+2,K)=TAUP_1(MM1+1,MM2+1,K)
         DO I=1,MM1+2
            TAUR_1(I,1,K)=TAUR_1(I,2,K)
            TAUP_1(I,1,K)=TAUP_1(I,2,K)
         ENDDO
         DO I=1,MM2+2
            TAUR_1(1,I,K)=TAUR_1(2,I,K)
            TAUP_1(1,I,K)=TAUP_1(2,I,K)
         ENDDO
         DO I=1,MM2+2
            TAUR_1(MM1+2,I,K)=TAUR_1(MM1+1,I,K)
            TAUP_1(MM1+2,I,K)=TAUP_1(MM1+1,I,K)
         ENDDO
         DO I=1,MM1+2
            TAUR_1(I,MM2+2,K)=TAUR_1(I,MM2+1,K)
            TAUP_1(I,MM2+2,K)=TAUP_1(I,MM2+1,K)
         ENDDO
      ENDDO

**********************************************************************
C     Integration of the ADT equations by 8th order Runge-Kutta method
********************************************************************** 

      include 'integration.h'

      STOP
      END

*************************************************SUBROUTINES**************************************


***********************************
C      Subroutine 'AMAT'
***********************************

**********************************************************************
*This subroutine returns numerically calculated ADT matrix (AA) for a
*particular set of ADT angles (Y).
**********************************************************************

      SUBROUTINE AMAT(Y,AA)
      IMPLICIT NONE
      INTEGER N,NN
      PARAMETER (N='NO OF STATE')
      PARAMETER (NN='NO OF TAU')
      INTEGER I,J,K,L,K1,I1
      REAL*8 ADD
      REAL*8 Y(NN),A(NN,N,N),AA1(N,N),AA(N,N)
 
      DO K=1,NN
        DO I=1,N
          DO J=1,N
           A(K,I,J)=0.0D0
          ENDDO
          A(K,I,I)=1.0D0
        ENDDO
      ENDDO
 
      K1=0
 
      DO J=2,N
        DO I=1,J-1

          K1=K1+1      !i.e.A12 is denoted as A1, A13 is denoted as A2, A23 is denoted as A3 etc.

          A(K1,I,I)=DCOS(Y(K1))
          A(K1,J,J)=DCOS(Y(K1))          
          A(K1,I,J)=DSIN(Y(K1))
          A(K1,J,I)=-DSIN(Y(K1))

        ENDDO
      ENDDO
 
      DO I=1,N
        DO J=1,N
         AA(I,J)=0.0D0
        ENDDO
         AA(I,I)=1.0D0
      ENDDO
 
      DO I1=1,NN
 
        DO I=1,N
           DO J=1,N
              ADD=0.0D0                    
              DO K=1,N
                ADD=ADD+AA(I,K)*A(I1,K,J)    
              ENDDO
              AA1(I,J)=ADD
           ENDDO
        ENDDO
 
        AA = AA1
 
      ENDDO
 
      RETURN
      END       

***********************************
C      Subroutine 'AMATSPL'
***********************************

***********************************************************************************************************
*This subroutine is implemented not only to generate the complete ADT matrix (AA), but also two sets of 
*partially multiplied ADT matrices (AFOR, ABACK). Adiabatic to diabatic transformation matrix
*is generated by multiplying elementary rotation matrices in a definite order, but partial ADT matrices 
*can be constructed by collecting one or more [2,3,4,.....,(N-1); N = number of elementary matrices] 
*matrices from that set and multiplying them in the same order as the parent one.
***********************************************************************************************************

      SUBROUTINE AMATSPL(Y,AA,AFOR,ABACK)
      IMPLICIT NONE
      INTEGER N,NN
      PARAMETER (N='NO OF STATE')
      PARAMETER (NN='NO OF TAU')
      INTEGER I,J,K,L,K1,I1
      REAL*8 ADD
      REAL*8 Y(NN),A(NN,N,N),AFOR(0:NN,N,N),ABACK(NN+1,N,N)
      REAL*8 AA1(N,N),AA2(N,N),AA(N,N)
 
      DO K=1,NN
        DO I=1,N
          DO J=1,N
           A(K,I,J)=0.0D0
          ENDDO
          A(K,I,I)=1.0D0
        ENDDO
      ENDDO
 
      K1=0
 
      DO J=2,N
        DO I=1,J-1

          K1=K1+1      !i.e.A12 is denoted as A1, A13 is denoted as A2, A23 is denoted as A3 etc.

          A(K1,I,I)=DCOS(Y(K1))
          A(K1,J,J)=DCOS(Y(K1))          
          A(K1,I,J)=DSIN(Y(K1))
          A(K1,J,I)=-DSIN(Y(K1))

        ENDDO
      ENDDO
 
      DO I=1,N
        DO J=1,N
         AA(I,J)=0.0D0
        ENDDO
         AA(I,I)=1.0D0
      ENDDO

      DO I1=1,NN
 
        DO I=1,N
           DO J=1,N
              ADD=0.0D0                    
              DO K=1,N
                ADD=ADD+AA(I,K)*A(I1,K,J)    
              ENDDO
              AA1(I,J)=ADD
           ENDDO
        ENDDO
 
        AA = AA1
        DO I=1,N
           DO J=1,N
              AFOR(I1,I,J)=AA1(I,J)
           ENDDO
        ENDDO
 
      ENDDO

      DO I=1,N
        DO J=1,N
         AA(I,J)=0.0D0
        ENDDO
         AA(I,I)=1.0D0
      ENDDO
 
      DO I1=NN,1,-1
 
        DO I=1,N
           DO J=1,N
              ADD=0.0D0                    
              DO K=1,N
                ADD=ADD+A(I1,I,K)*AA(K,J)    
              ENDDO
              AA2(I,J)=ADD
           ENDDO
        ENDDO
 
        AA = AA2
        DO I=1,N
           DO J=1,N
              ABACK(I1,I,J)=AA2(I,J)
           ENDDO
        ENDDO
 
      ENDDO

      DO I=1,N
         DO J=1,N
            AFOR(0,I,J)=0.0D0
         ENDDO
         AFOR(0,I,I)=1.0D0
      ENDDO

      DO I=1,N
         DO J=1,N
            ABACK(NN+1,I,J)=0.0D0
         ENDDO
         ABACK(NN+1,I,I)=1.0D0
      ENDDO

      RETURN
      END

************************************************************************************************************
C      Subroutine 'BCUCOF' (taken from 'W. H. Press, S. A. Teukolsky, W.T. Vetterling, B. P. Flannery, 
C      Numerical Recipes in FORTRAN, Cambridge University Press, A-229, DSIDC Industrial Park, Narela, 
C      Delhi-110040, 2000.')
************************************************************************************************************

******************************************************************************************************
*While integrating the differential equations by 8th order Runge-Kutta method, NACT values for the 
*intermediate geometries between two grid points are required and therefore, bi-cubic interpolation 
*is adopted to dig out the magnitude of NACTs at those unknown points. This subroutine is used as a 
*secondary subroutine in 'INTERPOL'.
******************************************************************************************************

      SUBROUTINE BCUCOF (Y,Y1,Y2,Y12,D1,D2,C)
      REAL*8 D1,D2,C(4,4),Y(4),Y1(4),Y12(4),Y2(4)
      INTEGER I,J,K,L
      REAL*8 D1D2,XX,CL(16),WT(16,16),X(16)
      SAVE WT
      DATA WT/1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4,8*0,3,0,-9,6,-2,0,6,-4
     * ,10*0,9,-6,2*0,-6,4,2*0,3,-2,6*0,-9,6,2*0,6,-4
     * ,4*0,1,0,-3,2,-2,0,6,-4,1,0,-3,2,8*0,-1,0,3,-2,1,0,-3,2
     * ,10*0,-3,2,2*0,3,-2,6*0,3,-2,2*0,-6,4,2*0,3,-2
     * ,0,1,-2,1,5*0,-3,6,-3,0,2,-4,2,9*0,3,-6,3,0,-2,4,-2
     * ,10*0,-3,3,2*0,2,-2,2*0,-1,1,6*0,3,-3,2*0,-2,2
     * ,5*0,1,-2,1,0,-2,4,-2,0,1,-2,1,9*0,-1,2,-1,0,1,-2,1
     * ,10*0,1,-1,2*0,-1,1,6*0,-1,1,2*0,2,-2,2*0,-1,1/
      D1D2=D1*D2
      DO I=1,4
      X(I)=Y(I)
      X(I+4)=Y1(I)*D1
      X(I+8)=Y2(I)*D2
      X(I+12)=Y12(I)*D1D2
      ENDDO
      DO I=1,16
      XX=0.0D0
      DO K=1,16
      XX=XX+WT(I,K)*X(K)
      ENDDO
      CL(I)=XX
      ENDDO
      L=0
      DO I=1,4
      DO J=1,4
      L=L+1
      C(I,J)=CL(L)
      ENDDO
      ENDDO
      RETURN
      END
   
************************************************************************************************************
C      Subroutine 'BCUINT' (taken from 'W. H. Press, S. A. Teukolsky, W.T. Vetterling, B. P. Flannery, 
C      Numerical Recipes in FORTRAN, Cambridge University Press, A-229, DSIDC Industrial Park, Narela, 
C      Delhi-110040, 2000.')
************************************************************************************************************

******************************************************************************************************
*While integrating the differential equations by 8th order Runge-Kutta method, NACT values for the 
*intermediate geometries between two grid points are required and therefore, bi-cubic interpolation 
*is adopted to dig out the magnitude of NACTs at those unknown points. This subroutine is used as a 
*secondary subroutine in 'INTERPOL'.
******************************************************************************************************

      SUBROUTINE BCUINT(Y,Y1,Y2,Y12,X1L,X1U,X2L,X2U,X1,X2,ANSY,
     * ANSY1,ANSY2)
      REAL*8 ANSY,ANSY1,ANSY2,X1,X1L,X1U,X2,X2L,X2U,Y(4),Y1(4),
     * Y12(4),Y2(4)
      INTEGER I
      REAL*8 T,U,C(4,4)
      CALL BCUCOF (Y,Y1,Y2,Y12,X1U-X1L,X2U-X2L,C)
      IF(X1U.EQ.X1L.OR.X2U.EQ.X2L)PAUSE 'bad input in bcuint'
      T=(X1-X1L)/(X1U-X1L)
      U=(X2-X2L)/(X2U-X2L)
      ANSY=0.0D0
      ANSY1=0.0D0
      ANSY2=0.0D0
      DO I=4,1,-1
      ANSY=T*ANSY+((C(I,4)*U+C(I,3))*U+C(I,2))*U+C(I,1)
      ANSY2=T*ANSY2+(3.0D0*C(I,4)*U+2.0D0*C(I,3))*U+C(I,2)
      ANSY1=U*ANSY1+(3.0D0*C(4,I)*T+2.0D0*C(3,I))*T+C(2,I)
      ENDDO
      ANSY1=ANSY1/(X1U-X1L)
      ANSY2=ANSY2/(X2U-X2L)
      RETURN
      END

***********************************
C      Subroutine 'FEXP'
***********************************

******************************************************************************************
*This subroutine returns the magnitude of $^{N}$C$_{2}$ numbers of ADT angles at a 
*specific nuclear geometry using a pre-calculated/predefined set of ADT angles (initial 
*value). This routine is implemented during integration along second coordinate when the 
*first coordinate is fixed at a definite value 
******************************************************************************************

      SUBROUTINE FEXP(C2,Y,FP)
      IMPLICIT NONE
      INTEGER N,MM1,MM2,I,J,NNN
      PARAMETER (MM1='C1 LIMIT',MM2='C2 LIMIT')
      PARAMETER (N='NO OF STATE',NNN='NO OF TAU')
      REAL*8 C1,C2
      REAL*8 DSEC
      REAL*8 RR,PH
      REAL*8 TAUP
      REAL*8 TAU(NNN)
      REAL*8 Y(NNN),FP(NNN)
      REAL*8 SOL(NNN)
      REAL*8 G(NNN,NNN),GI(NNN,NNN)
      REAL*8 A(NNN,N,N)
      COMMON/RHH/C1
      COMMON/GRID/RR(MM1+2),PH(MM2+2)
      COMMON/TPP/TAUP(MM1+2,MM2+2,NNN)

      CALL INTERPOL(MM1,MM2,NNN,RR,PH,TAUP,C1,C2,TAU)
      CALL RES(Y,TAU,FP)
 
      RETURN
      END

***********************************
C      Subroutine 'FEXR'
***********************************

******************************************************************************************
*This subroutine returns the magnitude of $^{N}$C$_{2}$ numbers of ADT angles at a 
*specific nuclear geometry using a pre-calculated/predefined set of ADT angles (initial 
*value). This routine is implemented during integration along first coordinate when the 
*second coordinate is fixed at a definite value 
******************************************************************************************

      SUBROUTINE FEXR(C1,Y,FR)
      IMPLICIT NONE
      INTEGER N,MM1,MM2,I,J,NNN
      PARAMETER (MM1='C1 LIMIT',MM2='C2 LIMIT')
      PARAMETER (N='NO OF STATE',NNN='NO OF TAU')
      REAL*8 C1,C2
      REAL*8 DSEC
      REAL*8 RR,PH
      REAL*8 TAUR
      REAL*8 TAU(NNN)
      REAL*8 Y(NNN),FR(NNN)
      REAL*8 SOL(NNN)
      REAL*8 G(NNN,NNN),GI(NNN,NNN)
      REAL*8 A(NNN,N,N)
      COMMON/PHH/C2
      COMMON/GRID/RR(MM1+2),PH(MM2+2)
      COMMON/TRR/TAUR(MM1+2,MM2+2,NNN)    
 
      CALL INTERPOL(MM1,MM2,NNN,RR,PH,TAUR,C1,C2,TAU)
      CALL RES(Y,TAU,FR)
 
      RETURN
      END
 
***********************************
C      Subroutine 'GRADCOMAT'
***********************************

***************************************************************************************
*This subprogram is implemented for calculating the elements of coefficient matrix 
*of gradient of ADT angles
***************************************************************************************

      SUBROUTINE GRADCOMAT(Y,AFOR,ABACK,G)
      IMPLICIT NONE
      INTEGER N,NN
      PARAMETER (N='NO OF STATE')
      PARAMETER (NN='NO OF TAU')
      INTEGER I,J,K,L,K1,I1,COUNTER
      REAL*8 ADD,START
      REAL*8 Y(NN)
      REAL*8 AFOR(0:NN,N,N),ABACK(NN+1,N,N),ADIFF(NN,N,N)
      REAL*8 AA1(N,N),AA(N,N),G(NN,NN)
 
      DO K=1,NN
        DO I=1,N
          DO J=1,N
            ADIFF(K,I,J)=0.0D0
          ENDDO
        ENDDO
      ENDDO
 
      K1=0
 
      DO J=2,N
        DO I=1,J-1

          K1=K1+1      !i.e.A12 is denoted as A1, A13 is denoted as A2, A23 is denoted as A3 etc.

          ADIFF(K1,I,I)=-DSIN(Y(K1))
          ADIFF(K1,J,J)=-DSIN(Y(K1))          
          ADIFF(K1,I,J)=DCOS(Y(K1))
          ADIFF(K1,J,I)=-DCOS(Y(K1))

        ENDDO
      ENDDO

      DO I1=1,NN

         DO I=1,N
            DO J=1,N
              ADD=0.0D0                    
              DO K=1,N
                   ADD=ADD+AFOR(I1-1,I,K)*ADIFF(I1,K,J)   
              ENDDO
              AA1(I,J)=ADD
            ENDDO
         ENDDO

         DO I=1,N
            DO J=1,N
              ADD=0.0D0                    
              DO K=1,N
                   ADD=ADD+AA1(I,K)*ABACK(I1+1,K,J)   
              ENDDO
              AA(I,J)=ADD
            ENDDO
         ENDDO

        COUNTER = 0
        DO I=2,N
          DO J=1,I-1
             COUNTER = COUNTER + 1
             G(COUNTER,I1) = AA(I,J)
          ENDDO
        ENDDO

      ENDDO
 
      RETURN
      END       

************************************************************************************************************
C      Subroutine 'INTERPOL' (taken from 'W. H. Press, S. A. Teukolsky, W.T. Vetterling, B. P. Flannery, 
C      Numerical Recipes in FORTRAN, Cambridge University Press, A-229, DSIDC Industrial Park, Narela, 
C      Delhi-110040, 2000.')
************************************************************************************************************

******************************************************************************************************
*While integrating the differential equations by 8th order Runge-Kutta method, NACT values for the 
*intermediate geometries between two grid points are required and therefore, bi-cubic interpolation 
*is adopted to dig out the magnitude of NACTs at those unknown points. This subroutine is used to  
*perform this interpolation.
******************************************************************************************************

      SUBROUTINE INTERPOL(MM1,MM2,NN,RR,PH,TAU,X1,Y1,TOUT)
      IMPLICIT NONE
      INTEGER MM1,MM2,NN
      INTEGER K
      REAL*8 RR(MM1+2),PH(MM2+2)
      REAL*8 TAU(MM1+2,MM2+2,NN)
      REAL*8 X1,Y1
      REAL*8 TOUT(NN)
      REAL*8 YF(4),YF1(4),YF2(4),YF12(4)
      REAL*8 DX,DY
      INTEGER II1,II2,JJ1,JJ2
      REAL*8 X1L,X1U,X2L,X2U
      REAL*8 ANSY,ANSY1,ANSY2
      DX=RR(2)-RR(1)
      DY=PH(2)-PH(1)
      CALL LOCATE(RR,MM1+2,X1,II1)
      II2=II1+1
      X1L=RR(II1)
      X1U=RR(II2)
      CALL LOCATE(PH,MM2+2,Y1,JJ1)
      JJ2=JJ1+1
      X2L=PH(JJ1)
      X2U=PH(JJ2)
      WRITE(8,55)II1,II2,JJ1,JJ2
      DO K=1,NN
C     ZEROTH DERIVATIVES
      YF(1)=TAU(II1,JJ1,K)
      YF(2)=TAU(II2,JJ1,K)
      YF(3)=TAU(II2,JJ2,K)
      YF(4)=TAU(II1,JJ2,K)
C     FIRST DERIVATIVES
      IF(II1.EQ.1)THEN
C     Forward Differencing(f.d)
      YF1(1)=(TAU(II1+1,JJ1,K)-TAU(II1,JJ1,K))/DX
      YF1(4)=(TAU(II1+1,JJ2,K)-TAU(II1,JJ2,K))/DX
      ELSE
C     Central Differencing(c.d)
      YF1(1)=0.5D0*(TAU(II1+1,JJ1,K)-TAU(II1-1,JJ1,K))/DX
      YF1(4)=0.5D0*(TAU(II1+1,JJ2,K)-TAU(II1-1,JJ2,K))/DX
      ENDIF
      IF(II2.EQ.MM1+2)THEN
C     Backward Differencing(b.d)
      YF1(2)=(TAU(II2,JJ1,K)-TAU(II2-1,JJ1,K))/DX
      YF1(3)=(TAU(II2,JJ2,K)-TAU(II2-1,JJ2,K))/DX
      ELSE
C     c.d
      YF1(2)=0.5D0*(TAU(II2+1,JJ1,K)-TAU(II2-1,JJ1,K))/DX
      YF1(3)=0.5D0*(TAU(II2+1,JJ2,K)-TAU(II2-1,JJ2,K))/DX
      ENDIF
      IF(JJ1.EQ.1)THEN
C     f.d
      YF2(1)=(TAU(II1,JJ1+1,K)-TAU(II1,JJ1,K))/DY
      YF2(2)=(TAU(II2,JJ1+1,K)-TAU(II2,JJ1,K))/DY
      ELSE
C     c.d
      YF2(1)=0.5D0*(TAU(II1,JJ1+1,K)-TAU(II1,JJ1-1,K))/DY
      YF2(2)=0.5D0*(TAU(II2,JJ1+1,K)-TAU(II2,JJ1-1,K))/DY
      ENDIF
      IF(JJ2.EQ.MM2+2)THEN
C     b.d
      YF2(3)=(TAU(II2,JJ2,K)-TAU(II2,JJ2-1,K))/DY
      YF2(4)=(TAU(II1,JJ2,K)-TAU(II1,JJ2-1,K))/DY
      ELSE
C     c.d
      YF2(3)=0.5D0*(TAU(II2,JJ2+1,K)-TAU(II2,JJ2-1,K))/DY
      YF2(4)=0.5D0*(TAU(II1,JJ2+1,K)-TAU(II1,JJ2-1,K))/DY
      ENDIF
C     SECOND DERIVATIVES
      IF(II1.EQ.1.AND.JJ1.EQ.1)THEN      
C     f.d & f.d
      YF12(1)=((TAU(II1+1,JJ1+1,K)-TAU(II1,JJ1+1,K))
     1-(TAU(II1+1,JJ1,K)-TAU(II1,JJ1,K)))/DX/DY
      ELSEIF(II1.EQ.1.AND.JJ1.NE.1)THEN      
C     f.d & c.d
      YF12(1)=0.5D0*((TAU(II1+1,JJ1+1,K)-TAU(II1,JJ1+1,K))
     1-(TAU(II1+1,JJ1-1,K)-TAU(II1,JJ1-1,K)))/DX/DY
      ELSEIF(II1.NE.1.AND.JJ1.EQ.1)THEN      
C     c.d & f.d
      YF12(1)=0.5D0*((TAU(II1+1,JJ1+1,K)-TAU(II1-1,JJ1+1,K))
     1-(TAU(II1+1,JJ1,K)-TAU(II1-1,JJ1,K)))/DX/DY
      ELSE
C     c.d & c.d
      YF12(1)=0.25D0*((TAU(II1+1,JJ1+1,K)-TAU(II1-1,JJ1+1,K))
     1-(TAU(II1+1,JJ1-1,K)-TAU(II1-1,JJ1-1,K)))/DX/DY
      ENDIF
      IF(II2.EQ.MM1+2.AND.JJ1.EQ.1)THEN      
C     b.d & f.d
      YF12(2)=((TAU(II2,JJ1+1,K)-TAU(II2-1,JJ1+1,K))
     1-(TAU(II2,JJ1,K)-TAU(II2-1,JJ1,K)))/DX/DY
      ELSEIF(II2.EQ.MM1+2.AND.JJ1.NE.1)THEN      
C     b.d & c.d
      YF12(2)=0.5D0*((TAU(II2,JJ1+1,K)-TAU(II2-1,JJ1+1,K))
     1-(TAU(II2,JJ1-1,K)-TAU(II2-1,JJ1-1,K)))/DX/DY
      ELSEIF(II2.NE.MM1+2.AND.JJ1.EQ.1)THEN      
C     c.d & f.d
      YF12(2)=0.5D0*((TAU(II2+1,JJ1+1,K)-TAU(II2-1,JJ1+1,K))
     1-(TAU(II2+1,JJ1,K)-TAU(II2-1,JJ1,K)))/DX/DY
      ELSE
C     c.d & c.d
      YF12(2)=0.25D0*((TAU(II2+1,JJ1+1,K)-TAU(II2-1,JJ1+1,K))
     1-(TAU(II2+1,JJ1-1,K)-TAU(II2-1,JJ1-1,K)))/DX/DY
      ENDIF
      IF(II2.EQ.MM1+2.AND.JJ2.EQ.MM2+2)THEN      
C     b.d & b.d
      YF12(3)=((TAU(II2,JJ2,K)+TAU(II2-1,JJ2,K))
     1-(TAU(II2,JJ2-1,K)-TAU(II2-1,JJ2-1,K)))/DX/DY
      ELSEIF(II2.EQ.MM1+2.AND.JJ2.NE.MM2+2)THEN      
C     b.d & c.d
      YF12(3)=0.5D0*((TAU(II2,JJ2+1,K)+TAU(II2-1,JJ2+1,K))
     1-(TAU(II2,JJ2-1,K)-TAU(II2-1,JJ2-1,K)))/DX/DY
      ELSEIF(II2.NE.MM1+2.AND.JJ2.EQ.MM2+2)THEN      
C     c.d & b.d
      YF12(3)=0.5D0*((TAU(II2+1,JJ2,K)+TAU(II2-1,JJ2,K))
     1-(TAU(II2+1,JJ2-1,K)-TAU(II2-1,JJ2-1,K)))/DX/DY
      ELSE
C     c.d & c.d
      YF12(3)=0.25D0*((TAU(II2+1,JJ2+1,K)+TAU(II2-1,JJ2+1,K))
     1-(TAU(II2+1,JJ2-1,K)-TAU(II2-1,JJ2-1,K)))/DX/DY
      ENDIF
C     f.d & b.d
      IF(II1.EQ.1.AND.JJ2.EQ.MM2+2)THEN      
      YF12(4)=((TAU(II1+1,JJ2,K)-TAU(II1,JJ2,K))
     1-(TAU(II1+1,JJ2-1,K)-TAU(II1,JJ2-1,K)))/DX/DY
      ELSEIF(II1.EQ.1.AND.JJ2.NE.MM2+2)THEN      
C     f.d & c.d
      YF12(4)=0.5D0*((TAU(II1+1,JJ2+1,K)-TAU(II1,JJ2+1,K))
     1-(TAU(II1+1,JJ2-1,K)-TAU(II1,JJ2-1,K)))/DX/DY
      ELSEIF(II1.NE.1.AND.JJ2.EQ.MM2+2)THEN      
C     c.d & b.d
      YF12(4)=0.5D0*((TAU(II1+1,JJ2,K)-TAU(II1-1,JJ2,K))
     1-(TAU(II1+1,JJ2-1,K)-TAU(II1-1,JJ2-1,K)))/DX/DY
      ELSE
C     c.d & c.d
      YF12(4)=0.25D0*((TAU(II1+1,JJ2+1,K)-TAU(II1-1,JJ2+1,K))
     1-(TAU(II1+1,JJ2-1,K)-TAU(II1-1,JJ2-1,K)))/DX/DY
      ENDIF
      CALL BCUINT(YF,YF1,YF2,YF12,X1L,X1U,X2L,X2U,X1,Y1,ANSY,
     *ANSY1,ANSY2)
      TOUT(K)=ANSY
      ENDDO
  55  FORMAT(1X,6I6)
  66  FORMAT(1X,10F12.6)
      RETURN
      END

***********************************
C      Subroutine 'INVERSE'
***********************************

*********************************************************************************
*The inverse of coefficient matrix of gradient of ADT angles is
*evaluated by Gauss-Jordon Method.
*********************************************************************************
      
      SUBROUTINE INVERSE(G,NN,GI)
      IMPLICIT NONE
      INTEGER NN
      INTEGER I,J,II,JJ,IRANK,IROW
      REAL*8 G(NN,NN),A(NN,NN)
      REAL*8 EN(NN,NN),GI(NN,NN)
      REAL*8 ZERO,TMP,FC

      ZERO = 1.0d-20         ! defines a numerical zero for the program

      A = G

      DO I = 1,NN
        DO J = 1,NN
          EN(I,J) = 0.0D0
          IF (I .EQ. J) EN(I,J) = 1.0D0
        ENDDO
      ENDDO
 
      GI = EN

      IRANK = NN             ! 'irank' is the dimension of the matrix to start with

C-------------------------------------------------------------
C     "Check whether the pivot element of the matrix is zero.
C-------------------------------------------------------------

 111  IROW = NN+1-IRANK

      IF (DABS(A(IROW,IROW)) .LE. ZERO) THEN          
         J = 0 
         DO I = IROW + 1,NN        
           IF (DABS(A(I,IROW)) .GT. ZERO) THEN 
              J = I
              EXIT
           ENDIF
         ENDDO

         IF (J .EQ. 0) THEN   ! this means no such row is found 
           WRITE (6,*) "Inverse of the given matrix does not exist!"
           GI = 0.0D0
           GOTO 222                   ! 222 --> return
         ELSE                         ! 'j'th row found for interchange
           DO I = 1,NN                
             TMP = A(IROW,I)
             A(IROW,I) = A(J,I)
             A(J,I) = TMP
             TMP = GI(IROW,I)
             GI(IROW,I) = GI(J,I)
             GI(J,I) = TMP
           ENDDO
         ENDIF

      ENDIF

      FC = A(IROW,IROW)

      DO I = 1,NN          
        A(IROW,I) = A(IROW,I)/FC
        GI(IROW,I) = GI(IROW,I)/FC    
      ENDDO

      DO I = 1,IROW - 1
        FC = A(I,IROW)             
        DO J = 1,NN         
          A(I,J) = A(I,J) - A(IROW,J)*FC
          GI(I,J) = GI(I,J) - GI(IROW,J)*FC
        ENDDO
      ENDDO

      DO I = IROW + 1,NN             
        FC = A(I,IROW)              
        DO J = 1,NN           
          A(I,J) = A(I,J) - A(IROW,J)*FC
          GI(I,J) = GI(I,J) - GI(IROW,J)*FC
        ENDDO
      ENDDO

      IF (IRANK .EQ. 1) THEN          
         GOTO 222
      ENDIF

      IRANK = IRANK - 1

      GOTO 111
        
 222  RETURN
      END

************************************************************************************************************
C      Subroutine 'LOCATE' (taken from 'W. H. Press, S. A. Teukolsky, W.T. Vetterling, B. P. Flannery, 
C      Numerical Recipes in FORTRAN, Cambridge University Press, A-229, DSIDC Industrial Park, Narela, 
C      Delhi-110040, 2000.')
************************************************************************************************************

******************************************************************************************************
*While integrating the differential equations by 8th order Runge-Kutta method, NACT values for the 
*intermediate geometries between two grid points are required and therefore, bi-cubic interpolation 
*is adopted to dig out the magnitude of NACTs at those unknown points. This subroutine is used as a 
*secondary subroutine in 'INTERPOL'.
******************************************************************************************************
      
      SUBROUTINE LOCATE(XX,N,X,J)
      IMPLICIT NONE
      INTEGER N
      INTEGER JL,JU,JM,J
      REAL*8 XX(N),X
      JL=0
      JU=N+1
 10   IF(JU-JL.GT.1)THEN
        JM=(JU+JL)/2
        IF((XX(N).GT.XX(1)).EQV.(X.GE.XX(JM)))THEN
           JL=JM
        ELSE
           JU=JM
        ENDIF
      GOTO 10
      ENDIF
      J=JL
      RETURN
      END

***********************************
C      Subroutine 'NEGTAU'
***********************************

******************************************************************************************
*This subprogram accepts the NACM as input and returns the negative form
*of the matrix at each and every grid point in nuclear CS.         
******************************************************************************************

      SUBROUTINE NEGTAU(TAU,TAUMAT)
      IMPLICIT NONE
      INTEGER N,NN
      PARAMETER (N='NO OF STATE')
      PARAMETER (NN='NO OF TAU')
      INTEGER I,J,COUNTER
      REAL*8 TAU(NN),TAUMAT(N,N)
 
      DO I=1,N
         TAUMAT(I,I) = 0.0D0
      ENDDO

      COUNTER = 0
      DO I=2,N
         DO J=1,I-1
              COUNTER = COUNTER + 1
              TAUMAT(I,J)=TAU(COUNTER)
         ENDDO
      ENDDO

      DO I=1,N-1
         DO J=I+1,N   
               TAUMAT(I,J)=-TAUMAT(J,I)
         ENDDO
      ENDDO

      RETURN

      END

***********************************
C      Subroutine 'RES'
***********************************

***************************************************************************************
*Magnitudes of gradient of ADT angles at every grid point in the nuclear
*CS is returned by this subroutine
***************************************************************************************

      SUBROUTINE RES(ANG,NACT,F)
      IMPLICIT NONE
      INTEGER N,NN,I,J,COUNTER
      PARAMETER (N='NO OF STATE')
      PARAMETER (NN='NO OF TAU') 
      REAL*8 ANG(NN),NACT(NN),VAL(NN),F(NN)
      REAL*8 G(NN,NN),GI(NN,NN) 
      REAL*8 AMATRIX(N,N),TMATRIX(N,N),PROD(N,N)
      REAL*8 A(NN,N,N),AFOR(0:NN,N,N),ABACK(NN+1,N,N)
      
      CALL AMATSPL(ANG,AMATRIX,AFOR,ABACK)
      CALL GRADCOMAT(ANG,AFOR,ABACK,G)
      CALL INVERSE(G,NN,GI)
      CALL NEGTAU(NACT,TMATRIX)

      PROD = MATMUL(TMATRIX,AMATRIX)
 
      COUNTER = 0
      DO I=2,N
        DO J=1,I-1
           COUNTER = COUNTER + 1
           VAL(COUNTER) = PROD(I,J)
        ENDDO
      ENDDO
      
      F = MATMUL(GI,VAL)

      RETURN
      END
        
***********************************
C      Subroutine 'WMAT'
***********************************

**********************************************************************************************
*The adiabatic potential energy (UMAT) is converted to the diabatic potential energy 
*matrix (WA) by performing similarity transformation with ADT matrix (AA) for the 
*user-defined number of coupled electronic states (NST).
**********************************************************************************************

      SUBROUTINE WMAT(NST,UMAT,AA,WA)
      IMPLICIT NONE 
      INTEGER NST
      INTEGER II,JJ,LL
      REAL*8 UMAT(NST,NST),AA(NST,NST)
      REAL*8 AADAG(NST,NST),BB(NST,NST)
      REAL*8 WA(NST,NST)
 
      DO II=1,NST
         DO JJ=1,NST
            AADAG(II,JJ)=AA(JJ,II)
         ENDDO
      ENDDO
 
      BB = MATMUL(UMAT,AA)
      WA = MATMUL(AADAG,BB)
 
      RETURN
      END

*************************************************FUNCTION**************************************

C     This function calculates secant of an angle

      FUNCTION DSEC(X)
      IMPLICIT NONE
      REAL*8 X,DSEC,F
      F=1/DCOS(X)
      DSEC=F
      RETURN
      END

***************************************************************************END******************************************************************
