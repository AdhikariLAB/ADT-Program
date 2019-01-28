***,NO2 PES and non-adiabatic coupling terms
memory,100,m
file,2,no2-q4.wfu,new;

basis=6-31G**;

ANG=0.529177;
HCROSS=0.063508;
NMASS=14.006700;
OMASS=15.999400;
CMINVTINV=0.001883651;

OMEGA1=848.07*CMINVTINV;
OMEGA2=1642.86*CMINVTINV;
OMEGA3=1959.73*CMINVTINV;

X10=0.000000;
Y10=0.000000;
Z10=-0.560164316;
X20=0.00000;
Y20=2.031159929;
Z20=0.245198367;
X30=0.00000;
Y30=-2.031159929;
Z30=0.245198367;

Q3=0;

RHOMIN=0.0;
RHOMAX=5.0;
DRHO=(RHOMAX-RHOMIN)/50.0;

PHIMIN=1.0 DEGREE;
PHIMAX=359.0 DEGREE;
DPHI=(PHIMAX-PHIMIN)/179.0;

PROC CALCULATE_COORD

x1=SQRT(HCROSS/NMASS)*(0.0000*Q1/SQRT(OMEGA1)+0.0000000*Q2/SQRT(OMEGA2)+0.0000000*Q3/SQRT(OMEGA3));


y1=SQRT(HCROSS/NMASS)*(0.0000000*Q1/SQRT(OMEGA1)+0.0000000*Q2/SQRT(OMEGA2)+0.8147216*Q3/SQRT(OMEGA3));


z1=SQRT(HCROSS/NMASS)*(0.6320366*Q1/SQRT(OMEGA1)+0.5441256*Q2/SQRT(OMEGA2)+0.0000000*Q3/SQRT(OMEGA3));


x2=SQRT(HCROSS/OMASS)*(0.0000*Q1/SQRT(OMEGA1)+0.0000000*Q2/SQRT(OMEGA2)+0.0000000*Q3/SQRT(OMEGA3));


y2=SQRT(HCROSS/OMASS)*(0.4613412*Q1/SQRT(OMEGA1)-0.5358772*Q2/SQRT(OMEGA2)-0.3811495*Q3/SQRT(OMEGA3));


z2=SQRT(HCROSS/OMASS)*(-0.2956843*Q1/SQRT(OMEGA1)-0.2545571*Q2/SQRT(OMEGA2)-0.1511272*Q3/SQRT(OMEGA3));


x3=SQRT(HCROSS/OMASS)*(-0.0000*Q1/SQRT(OMEGA1)+0.0000000*Q2/SQRT(OMEGA2)-0.0000*Q3/SQRT(OMEGA3));


y3=SQRT(HCROSS/OMASS)*(-0.463412*Q1/SQRT(OMEGA1)+0.5358772*Q2/SQRT(OMEGA2)-0.3811495*Q3/SQRT(OMEGA3));


z3=SQRT(HCROSS/OMASS)*(-0.2956843*Q1/SQRT(OMEGA1)-0.2545571*Q2/SQRT(OMEGA2)+0.1511272*Q3/SQRT(OMEGA3));

ENDPROC

!---------------------------------------
!EQUILIBRIUM POINT
!--------------------------------------

PHI=0.0
RHO=0.0

Q1=RHO*COS(PHI);
Q2=RHO*SIN(PHI);

CALCULATE_COORD

symmetry,nosym

geomtyp=xyz
geometry= {
3
This is the geometry input for NO2 in XYZ format
n,x1+x10,y1+y10,z1+z10
o,x2+x20,y2+y20,z2+z20
o,x3+x30,y3+y30,z3+z30
}

{uhf;accu,5;orbital,2200.2}
{multi;occ,13;closed,5; wf,23,1,1,0;state,2; maxiter,40;orbital,2140.2;}

show, energy
table, energy
save,no3-enrq3-q4.res
title 'results: Geom --- Equilibrium
{table,____; noprint,heading}

!----------------------
!Equilibrium ends here
!----------------------

!------------------------
!Normal grids start here
!-----------------------

do I=1,2

PHI=PHIMIN+(I-1)*DPHI
do J=1,10
RHO=RHOMIN+J*DRHO;

Q1=RHO*COS(PHI);
Q2=RHO*SIN(PHI);

CALCULATE_COORD

symmetry,nosym

geomtyp=xyz
geometry= {
3
This is the geometry input for NO2 in XYZ format
n,x1+x10,y1+y10,z1+z10
o,x2+x20,y2+y20,z2+z20
o,x3+x30,y3+y30,z3+z30

}


{multi;occ,13;closed,5; wf,23,1,1,0;state,2; maxiter,40;start,2140.2;orbital,2140.2;}


show, energy
table, energy
save,no3-enrq3-q4.res
title 'results: Geom --- RHO= $RHO PHI= $PHI
{table,____; noprint,heading}



{multi;occ,13;closed,5; wf,23,1,1,0;state,2;start,2140.2;orbital,2350.2;maxiter,40;
cpmcscf,nacm,1.1,2.1,accu=1.d-10,record=5101.1;}


show, energy
table, energy

force;nacm,5101.1;varsav
table,gradx,grady,gradz;
save,no3-analytic-nac.res;
title 'results: Tau_12 ---Geom --- RHO= $RHO PHI= $PHI
{table,____;  noprint,heading}

enddo
enddo

---
