***,NO3 PES and non-adiabatic coupling terms
file,2,no3_pes_tauq3-q4.wfu,new;

basis=6-31G**;

ANG=0.529177;
HCROSS=0.063508;
NMASS=14.006700;
OMASS=15.999400;
CMINVTINV=0.001883651;

OMEGA1=1047*CMINVTINV;
OMEGA2=745*CMINVTINV;
OMEGA3=1185*CMINVTINV;
OMEGA4=1185*CMINVTINV;
OMEGA5=476*CMINVTINV;
OMEGA6=476*CMINVTINV;

X10=0.000000*ANG;
Y10=0.000000*ANG;
Z10=0.000000*ANG;
X20=2.337215*ANG;
Y20=0.000000*ANG;
Z20=0.000000*ANG;
X30=-1.168607*ANG;
Y30=-2.024086*ANG;
Z30=0.000000*ANG;
X40=-1.168607*ANG;
Y40=2.024086*ANG;
Z40=0.000000*ANG;

Q1=0.00; 
Q2=0.00; 
Q5=0.00;
Q6=0.00;

RHOMIN=0.0;
RHOMAX=5.0;
DRHO=(RHOMAX-RHOMIN)/50.0;

PHIMIN=1.0 DEGREE;
PHIMAX=359.0 DEGREE;
DPHI=(PHIMAX-PHIMIN)/179.0;

PROC CALCULATE_COORD

ax1=SQRT(HCROSS/NMASS)*(0.0000102*Q1/SQRT(OMEGA1)+0.0000000*Q2/SQRT(OMEGA2)+0.0000000*Q3/SQRT(OMEGA3));
bx1=SQRT(HCROSS/NMASS)*(+0.7449584*Q4/SQRT(OMEGA4)+0.4681243*Q5/SQRT(OMEGA5)+0.0000000*Q6/SQRT(OMEGA6));
x1=ax1+bx1;

ay1=SQRT(HCROSS/NMASS)*(0.0000000*Q1/SQRT(OMEGA1)+0.0000000*Q2/SQRT(OMEGA2)+0.7448872*Q3/SQRT(OMEGA3));
by1=SQRT(HCROSS/NMASS)*(0.0000000*Q4/SQRT(OMEGA4)+0.0000000*Q5/SQRT(OMEGA5)-0.4682375*Q6/SQRT(OMEGA6));
y1=ay1+by1;

az1=SQRT(HCROSS/NMASS)*(0.0000000*Q1/SQRT(OMEGA1)+0.8798314*Q2/SQRT(OMEGA2)+0.0000000*Q3/SQRT(OMEGA3));
bz1=SQRT(HCROSS/NMASS)*(0.0000000*Q4/SQRT(OMEGA4)+0.0000000*Q5/SQRT(OMEGA5)+0.0000000*Q6/SQRT(OMEGA6));
z1=az1+bz1;

ax2=SQRT(HCROSS/OMASS)*(0.5774593*Q1/SQRT(OMEGA1)+0.0000000*Q2/SQRT(OMEGA2)+0.0000000*Q3/SQRT(OMEGA3));
bx2=SQRT(HCROSS/OMASS)*(-0.5394732*Q4/SQRT(OMEGA4)+0.3427462*Q5/SQRT(OMEGA5)+0.0000000*Q6/SQRT(OMEGA6));
x2=ax2+bx2;

ay2=SQRT(HCROSS/OMASS)*(0.0000000*Q1/SQRT(OMEGA1)+0.0000000*Q2/SQRT(OMEGA2)+0.0749380*Q3/SQRT(OMEGA3));
by2=SQRT(HCROSS/OMASS)*(0.0000000*Q4/SQRT(OMEGA4)+0.0000000*Q5/SQRT(OMEGA5)+0.6348331*Q6/SQRT(OMEGA6));
y2=ay2+by2;

az2=SQRT(HCROSS/OMASS)*(0.0000000*Q1/SQRT(OMEGA1)-0.2744083*Q2/SQRT(OMEGA2)+0.0000000*Q3/SQRT(OMEGA3));
bz2=SQRT(HCROSS/OMASS)*(0.0000000*Q4/SQRT(OMEGA4)+0.0000000*Q5/SQRT(OMEGA5)+0.0000000*Q6/SQRT(OMEGA6));
z2=az2+bz2;

ax3=SQRT(HCROSS/OMASS)*(-0.2887344*Q1/SQRT(OMEGA1)+0.0000000*Q2/SQRT(OMEGA2)-0.2660973*Q3/SQRT(OMEGA3));
bx3=SQRT(HCROSS/OMASS)*(-0.0787757*Q4/SQRT(OMEGA4)-0.3903747*Q5/SQRT(OMEGA5)-0.4233146*Q6/SQRT(OMEGA6));
x3=ax3+bx3;

ay3=SQRT(HCROSS/OMASS)*(-0.4999028*Q1/SQRT(OMEGA1)+0.0000000*Q2/SQRT(OMEGA2)-0.3859480*Q3/SQRT(OMEGA3));
by3=SQRT(HCROSS/OMASS)*(-0.2660775*Q4/SQRT(OMEGA4)+0.4234383*Q5/SQRT(OMEGA5)-0.0983620*Q6/SQRT(OMEGA6));
y3=ay3+by3;

az3=SQRT(HCROSS/OMASS)*(0.0000000*Q1/SQRT(OMEGA1)-0.2744083*Q2/SQRT(OMEGA2)+0.0000000*Q3/SQRT(OMEGA3));
bz3=SQRT(HCROSS/OMASS)*(0.0000000*Q4/SQRT(OMEGA4)+0.0000000*Q5/SQRT(OMEGA5)+0.0000000*Q6/SQRT(OMEGA6));
z3=az3+bz3;

ax4=SQRT(HCROSS/OMASS)*(-0.2887344*Q1/SQRT(OMEGA1)+0.0000000*Q2/SQRT(OMEGA2)+0.2660973*Q3/SQRT(OMEGA3));
bx4=SQRT(HCROSS/OMASS)*(-0.0787757*Q4/SQRT(OMEGA4)-0.3903747*Q5/SQRT(OMEGA5)+0.4233146*Q6/SQRT(OMEGA6));
x4=ax4+bx4;

ay4=SQRT(HCROSS/OMASS)*(0.4999028*Q1/SQRT(OMEGA1)+0.0000000*Q2/SQRT(OMEGA2)-0.3859480*Q3/SQRT(OMEGA3));
by4=SQRT(HCROSS/OMASS)*(0.2660775*Q4/SQRT(OMEGA4)-0.4234383*Q5/SQRT(OMEGA5)-0.0983620*Q6/SQRT(OMEGA6));
y4=ay4+by4;

az4=SQRT(HCROSS/OMASS)*(0.0000000*Q1/SQRT(OMEGA1)-0.2744083*Q2/SQRT(OMEGA2)+0.0000000*Q3/SQRT(OMEGA3));
bz4=SQRT(HCROSS/OMASS)*(0.0000000*Q4/SQRT(OMEGA4)+0.0000000*Q5/SQRT(OMEGA5)+0.0000000*Q6/SQRT(OMEGA6));
z4=az4+bz4;

ENDPROC

!---------------------------------------
!EQUILIBRIUM POINT
!--------------------------------------

PHI=0.0
RHO=0.0

Q3=RHO*COS(PHI);
Q4=RHO*SIN(PHI);

CALCULATE_COORD

symmetry,nosym

geomtyp=xyz
geometry= {
4
This is the geometry input for NO3 in XYZ format
n,x1+x10,y1+y10,z1+z10
o,x2+x20,y2+y20,z2+z20
o,x3+x30,y3+y30,z3+z30
o,x4+x40,y4+y40,z4+z40
}

{uhf;accu,5;orbital,2200.2}
{multi;occ,19;closed,11;wf,31,1,1,0;state,5;
       start,2200.2;orbital,2140.2;maxiter,40;}

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

DO I=1,60
PHI=PHIMIN+(I-1)*DPHI
!IF(I.EQ.25)THEN
!K=38
!ELSE
!K=1
!ENDIF
DO J=1,50
RHO=RHOMIN+J*DRHO;

Q3=RHO*COS(PHI);
Q4=RHO*SIN(PHI);

CALCULATE_COORD

symmetry,nosym

geomtyp=xyz
geometry= {
4
This is the geometry input for NO3 in XYZ format
n,x1+x10,y1+y10,z1+z10
o,x2+x20,y2+y20,z2+z20
o,x3+x30,y3+y30,z3+z30
o,x4+x40,y4+y40,z4+z40
}

if(j.eq.1)then
{multi;occ,19;closed,11;wf,31,1,1,0;state,5;
       start,2140.2;orbital,2140.2+j;maxiter,40;}
!elseif(i.eq.1.and.j.eq.k)then
!{multi;occ,19;closed,11;wf,31,1,1,0;state,5;
!       start,2140.2+(j-2);orbital,2140.2+j;maxiter,40;}
else
{multi;occ,19;closed,11;wf,31,1,1,0;state,5;
       start,2140.2+(j-1);orbital,2140.2+j;maxiter,40;}
endif

show, energy
table, energy
save,no3-enrq3-q4.res
title 'results: Geom --- RHO= $RHO PHI= $PHI
{table,____; noprint,heading}

{multi;occ,19;closed,11;wf,31,1,1,0;state,5;start,2140.2+j;orbital,2350.2;maxiter,40;
cpmcscf,nacm,1.1,2.1,accu=1.d-10,record=5101.1;
cpmcscf,nacm,1.1,3.1,accu=1.d-10,record=5102.1;
cpmcscf,nacm,1.1,4.1,accu=1.d-10,record=5103.1;
cpmcscf,nacm,1.1,5.1,accu=1.d-10,record=5104.1;
}

{multi;occ,19;closed,11;wf,31,1,1,0;state,5;start,2140.2+j;orbital,2350.2;maxiter,40;
cpmcscf,nacm,2.1,3.1,accu=1.d-10,record=5105.1;
cpmcscf,nacm,2.1,4.1,accu=1.d-10,record=5106.1;
cpmcscf,nacm,2.1,5.1,accu=1.d-10,record=5107.1;
}

{multi;occ,19;closed,11;wf,31,1,1,0;state,5;start,2140.2+j;orbital,2350.2;maxiter,40;
cpmcscf,nacm,3.1,4.1,accu=1.d-10,record=5108.1;
cpmcscf,nacm,3.1,5.1,accu=1.d-10,record=5109.1;
cpmcscf,nacm,4.1,5.1,accu=1.d-10,record=5110.1;
}

show, energy
table, energy

force;nacm,5101.1;varsav
table,gradx,grady,gradz;
save,no3-analytic-nac.res;
title 'results: Tau_12 ---Geom --- RHO= $RHO PHI= $PHI
{table,____;  noprint,heading}

force;nacm,5102.1;varsav
table,gradx,grady,gradz;
save,no3-analytic-nac.res;
title 'results: Tau_13 ---Geom --- RHO= $RHO PHI= $PHI
{table,____; noprint,heading}

force;nacm,5103.1;varsav
table,gradx,grady,gradz;
save,no3-analytic-nac.res;
title 'results: Tau_14 ---Geom --- RHO= $RHO PHI= $PHI
{table,____; noprint,heading}

 force;nacm,5104.1;varsav
 table,gradx,grady,gradz;
 save,no3-analytic-nac.res;
 title 'results: Tau_15 ---Geom --- RHO= $RHO PHI= $PHI
 {table,____; noprint,heading}

 force;nacm,5105.1;varsav
 table,gradx,grady,gradz;
 save,no3-analytic-nac.res;
 title 'results: Tau_23 ---Geom --- RHO= $RHO PHI= $PHI
 {table,____; noprint,heading}

 force;nacm,5106.1;varsav
 table,gradx,grady,gradz;
 save,no3-analytic-nac.res;
 title 'results: Tau_24 ---Geom --- RHO= $RHO PHI= $PHI
 {table,____; noprint,heading}

 force;nacm,5107.1;varsav
 table,gradx,grady,gradz;
 save,no3-analytic-nac.res;
 title 'results: Tau_25 ---Geom --- RHO= $RHO PHI= $PHI
 {table,____; noprint,heading}

 force;nacm,5108.1;varsav
 table,gradx,grady,gradz;
 save,no3-analytic-nac.res;
 title 'results: Tau_34 ---Geom --- RHO= $RHO PHI= $PHI
 {table,____; noprint,heading}

 force;nacm,5109.1;varsav
 table,gradx,grady,gradz;
 save,no3-analytic-nac.res;
 title 'results: Tau_35 ---Geom --- RHO= $RHO PHI= $PHI
 {table,____; noprint,heading}

 force;nacm,5110.1;varsav
 table,gradx,grady,gradz;
 save,no3-analytic-nac.res;
 title 'results: Tau_45 ---Geom --- RHO= $RHO PHI= $PHI
 {table,____; noprint,heading}

ENDDO
ENDDO

!------------------------
!Normal grids start here
!-----------------------
 ---
