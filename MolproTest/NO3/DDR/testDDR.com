***,NO3 PES and non-adiabatic coupling terms
memory,500,m
file,2,no3_pes-tau_q3-q4-new.wfu;
gthresh,twoint=1.d-15,oneint=1.d-15,energy=1.d-10,gradient=1.d-6,prefac=1.d-15

basis=cc-pVDZ;

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

DR=0.001;
DP=0.01;

Q1=0.00; 
Q2=0.00; 
Q5=0.00;
Q6=0.00;

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

!------------------------------------------------------------
!EQUILIBRIUM POINT
!------------------------------------------------------------

!Q3=0.0
!Q4=0.0

!CALCULATE_COORD

!symmetry,nosym

!geomtyp=xyz
!geometry= {
!4
!This is the geometry input for NO3 in XYZ format
!n,x1+x10,y1+y10,z1+z10
!o,x2+x20,y2+y20,z2+z20
!o,x3+x30,y3+y30,z3+z30
!o,x4+x40,y4+y40,z4+z40
!}

!{uhf;accu,5;orbital,2200.2}
!{multi;occ,19;closed,11;wf,31,1,1,0;state,5;
!       start,2200.2;orbital,2140.2;maxiter,40;}

!show, energy
!table, energy

!{mrci;occ,19;core,11;closed,11;wf,31,1,1,0;state,5}

!show, energy
!table, energy
!save,no3-enrq3-q4.res
!title 'results: EQUILIBRIUM'
!{table,____; noprint,heading}

!------------------------------------------------------------
!GRIDS
!------------------------------------------------------------

PHIMIN=0.0
PHIMAX=120.0
DPHI=2.0
RHOMIN=0.0
RHOMAX=5.0
DRHO=0.1


PHI=1
RHO=0.1

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


{multi;occ,19;closed,11;wf,31,1,1,0;state,5;
       start,2140.2;orbital,2140.2;maxiter,40;}

show, energy
table, energy

{mrci;occ,19;core,11;closed,11;wf,31,1,1,0;state,5;MAXITER,99,99;
             save,6000.2;dm,8000.2}

show, energy
table, energy
save,no3-enrq3-q4.res
title 'results: Geom --- RHO= $RHO PHI= $PHI
{table,____; noprint,heading}

!Calculate TAUR

RHO11=RHO+dr
Q3=RHO11*COS(PHI);
Q4=RHO11*SIN(PHI);

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


{multi;occ,19;closed,11;wf,31,1,1,0;state,5;maxiter,40;start,2140.2;
          orbital,2241.2;}

{mrci;occ,19;core,11;closed,11;wf,31,1,1,0;state,5;MAXITER,99,99;save,6001.2}

{ci;trans,6000.2,6001.2;
dm,8100.2}

RHO11=RHO-dr
Q3=RHO11*COS(PHI);
Q4=RHO11*SIN(PHI);

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

{multi;occ,19;closed,11;wf,31,1,1,0;state,5;maxiter,40;start,2140.2;
          orbital,2242.2;}

{mrci;occ,19;core,11;closed,11;wf,31,1,1,0;state,5;MAXITER,99,99;save,6002.2}

{ci;trans,6000.2,6002.2;
dm,8200.2}

{ddr,dr,2140.2,2241.2,8100.2;
state, 2.1,1.1}
nacme1p=nacme

{ddr,-dr,2140.2,2242.2,8200.2;
state, 2.1,1.1}
nacme1m=nacme

{ddr,dr,2140.2,2241.2,8100.2;
state, 3.1,1.1}
nacme2p=nacme

{ddr,-dr,2140.2,2242.2,8200.2;
state, 3.1,1.1}
nacme2m=nacme

{ddr,dr,2140.2,2241.2,8100.2;
state, 3.1,2.1}
nacme3p=nacme

{ddr,-dr,2140.2,2242.2,8200.2;
state, 3.1,2.1}
nacme3m=nacme

{ddr,2*dr
orbital,2140.2,2241.2,2242.2;
density,8000.2,8100.2,8200.2;
state, 2.1,1.1}
nacme2_1=nacme

{ddr,2*dr
orbital,2140.2,2241.2,2242.2;
density,8000.2,8100.2,8200.2;
state, 3.1,1.1}
nacme2_2=nacme

{ddr,2*dr
orbital,2140.2,2241.2,2242.2;
density,8000.2,8100.2,8200.2;
state, 3.1,2.1}
nacme2_3=nacme

nacmeav1=(nacme1p+nacme1m)*0.5
nacmeav2=(nacme2p+nacme2m)*0.5
nacmeav3=(nacme3p+nacme3m)*0.5

table,RHO,PHI,nacmeav1,nacmeav2,nacmeav3,nacme2_1,nacme2_2,nacme2_3
title,Non-adiabatic couplings for NO3
save,no3-nac-rho-q3-q4_123.res

{ddr,dr,2140.2,2241.2,8100.2;
state, 4.1,1.1}
nacme4p=nacme

{ddr,-dr,2140.2,2242.2,8200.2;
state, 4.1,1.1}
nacme4m=nacme

{ddr,dr,2140.2,2241.2,8100.2;
state, 5.1,1.1}
nacme5p=nacme

{ddr,-dr,2140.2,2242.2,8200.2;
state, 5.1,1.1}
nacme5m=nacme

{ddr,dr,2140.2,2241.2,8100.2;
state, 5.1,4.1}
nacme6p=nacme

{ddr,-dr,2140.2,2242.2,8200.2;
state, 5.1,4.1}
nacme6m=nacme

{ddr,2*dr
orbital,2140.2,2241.2,2242.2;
density,8000.2,8100.2,8200.2;
state, 4.1,1.1}
nacme2_4=nacme

{ddr,2*dr
orbital,2140.2,2241.2,2242.2;
density,8000.2,8100.2,8200.2;
state, 5.1,1.1}
nacme2_5=nacme

{ddr,2*dr
orbital,2140.2,2241.2,2242.2;
density,8000.2,8100.2,8200.2;
state, 5.1,4.1}
nacme2_6=nacme

nacmeav4=(nacme4p+nacme4m)*0.5
nacmeav5=(nacme5p+nacme5m)*0.5
nacmeav6=(nacme6p+nacme6m)*0.5

table,RHO,PHI,nacmeav4,nacmeav5,nacmeav6,nacme2_4,nacme2_5,nacme2_6
title,Non-adiabatic couplings for NO3
save,no3-nac-rho-q3-q4_145.res

{ddr,dr,2140.2,2241.2,8100.2;
state, 4.1,2.1}
nacme7p=nacme

{ddr,-dr,2140.2,2242.2,8200.2;
state, 4.1,2.1}
nacme7m=nacme

{ddr,dr,2140.2,2241.2,8100.2;
state, 5.1,2.1}
nacme8p=nacme

{ddr,-dr,2140.2,2242.2,8200.2;
state, 5.1,2.1}
nacme8m=nacme

{ddr,dr,2140.2,2241.2,8100.2;
state, 4.1,3.1}
nacme9p=nacme

{ddr,-dr,2140.2,2242.2,8200.2;
state, 4.1,3.1}
nacme9m=nacme

{ddr,dr,2140.2,2241.2,8100.2;
state, 5.1,3.1}
nacme10p=nacme

{ddr,-dr,2140.2,2242.2,8200.2;
state, 5.1,3.1}
nacme10m=nacme

{ddr,2*dr
orbital,2140.2,2241.2,2242.2;
density,8000.2,8100.2,8200.2;
state, 4.1,2.1}
nacme2_7=nacme

{ddr,2*dr
orbital,2140.2,2241.2,2242.2;
density,8000.2,8100.2,8200.2;
state, 5.1,2.1}
nacme2_8=nacme

{ddr,2*dr
orbital,2140.2,2241.2,2242.2;
density,8000.2,8100.2,8200.2;
state, 4.1,3.1}
nacme2_9=nacme

{ddr,2*dr
orbital,2140.2,2241.2,2242.2;
density,8000.2,8100.2,8200.2;
state, 5.1,3.1}
nacme2_10=nacme

nacmeav7=(nacme7p+nacme7m)*0.5
nacmeav8=(nacme8p+nacme8m)*0.5
nacmeav9=(nacme9p+nacme9m)*0.5
nacmeav10=(nacme10p+nacme10m)*0.5

table,RHO,PHI,nacmeav7,nacmeav8,nacmeav9,nacmeav10,nacme2_7,nacme2_8,nacme2_9,nacme2_10
title,Non-adiabatic couplings for NO3
save,no3-nac-rho-q3-q4_2345.res

!Calculate TAUP

PHI11=PHI+dp
Q3=RHO*COS(PHI11);
Q4=RHO*SIN(PHI11);

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

{multi;occ,19;closed,11;wf,31,1,1,0;state,5;maxiter,40;start,2140.2;
          orbital,2341.2;}

{mrci;occ,19;core,11;closed,11;wf,31,1,1,0;state,5;MAXITER,99,99;save,6001.2}

{ci;trans,6000.2,6001.2;
dm,8100.2}

PHI11=PHI-dp
Q3=RHO*COS(PHI11);
Q4=RHO*SIN(PHI11);

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

{multi;occ,19;closed,11;wf,31,1,1,0;state,5;maxiter,40;start,2140.2;
          orbital,2342.2;}

{mrci;occ,19;core,11;closed,11;wf,31,1,1,0;state,5;MAXITER,99,99;save,6002.2;}

{ci;trans,6000.2,6002.2;
dm,8200.2}

{ddr,dp,2140.2,2341.2,8100.2;
state, 2.1,1.1}
nacme1p=nacme

{ddr,-dp,2140.2,2342.2,8200.2;
state, 2.1,1.1}
nacme1m=nacme

{ddr,dp,2140.2,2341.2,8100.2;
state, 3.1,1.1}
nacme2p=nacme

{ddr,-dp,2140.2,2342.2,8200.2;
state, 3.1,1.1}
nacme2m=nacme

{ddr,dp,2140.2,2341.2,8100.2;
state, 3.1,2.1}
nacme3p=nacme

{ddr,-dp,2140.2,2342.2,8200.2;
state, 3.1,2.1}
nacme3m=nacme

{ddr,2*dp
orbital,2140.2,2341.2,2342.2;
density,8000.2,8100.2,8200.2;
state, 2.1,1.1}
nacme2_1=nacme

{ddr,2*dp
orbital,2140.2,2341.2,2342.2;
density,8000.2,8100.2,8200.2;
state, 3.1,1.1}
nacme2_2=nacme

{ddr,2*dp
orbital,2140.2,2341.2,2342.2;
density,8000.2,8100.2,8200.2;
state, 3.1,2.1}
nacme2_3=nacme

nacmeav1=(nacme1p+nacme1m)*0.5
nacmeav2=(nacme2p+nacme2m)*0.5
nacmeav3=(nacme3p+nacme3m)*0.5

table,RHO,PHI,nacmeav1,nacmeav2,nacmeav3,nacme2_1,nacme2_2,nacme2_3
title,Non-adiabatic couplings for NO3
save,no3-nac-phi-q3-q4_123.res

{ddr,dp,2140.2,2341.2,8100.2;
state, 4.1,1.1}
nacme4p=nacme

{ddr,-dp,2140.2,2342.2,8200.2;
state, 4.1,1.1}
nacme4m=nacme

{ddr,dp,2140.2,2341.2,8100.2;
state, 5.1,1.1}
nacme5p=nacme

{ddr,-dp,2140.2,2342.2,8200.2;
state, 5.1,1.1}
nacme5m=nacme

{ddr,dp,2140.2,2341.2,8100.2;
state, 5.1,4.1}
nacme6p=nacme

{ddr,-dp,2140.2,2342.2,8200.2;
state, 5.1,4.1}
nacme6m=nacme

{ddr,2*dp
orbital,2140.2,2341.2,2342.2;
density,8000.2,8100.2,8200.2;
state, 4.1,1.1}
nacme2_4=nacme

{ddr,2*dp
orbital,2140.2,2341.2,2342.2;
density,8000.2,8100.2,8200.2;
state, 5.1,1.1}
nacme2_5=nacme

{ddr,2*dp
orbital,2140.2,2341.2,2342.2;
density,8000.2,8100.2,8200.2;
state, 5.1,4.1}
nacme2_6=nacme

nacmeav4=(nacme4p+nacme4m)*0.5
nacmeav5=(nacme5p+nacme5m)*0.5
nacmeav6=(nacme6p+nacme6m)*0.5

table,RHO,PHI,nacmeav4,nacmeav5,nacmeav6,nacme2_4,nacme2_5,nacme2_6
title,Non-adiabatic couplings for NO3
save,no3-nac-phi-q3-q4_145.res

{ddr,dp,2140.2,2341.2,8100.2;
state, 4.1,2.1}
nacme7p=nacme

{ddr,-dp,2140.2,2342.2,8200.2;
state, 4.1,2.1}
nacme7m=nacme

{ddr,dp,2140.2,2341.2,8100.2;
state, 5.1,2.1}
nacme8p=nacme

{ddr,-dp,2140.2,2342.2,8200.2;
state, 5.1,2.1}
nacme8m=nacme

{ddr,dp,2140.2,2341.2,8100.2;
state, 4.1,3.1}
nacme9p=nacme

{ddr,-dp,2140.2,2342.2,8200.2;
state, 4.1,3.1}
nacme9m=nacme

{ddr,dp,2140.2,2341.2,8100.2;
state, 5.1,3.1}
nacme10p=nacme

{ddr,-dp,2140.2,2342.2,8200.2;
state, 5.1,3.1}
nacme10m=nacme

{ddr,2*dp
orbital,2140.2,2341.2,2342.2;
density,8000.2,8100.2,8200.2;
state, 4.1,2.1}
nacme2_7=nacme

{ddr,2*dp
orbital,2140.2,2341.2,2342.2;
density,8000.2,8100.2,8200.2;
state, 5.1,2.1}
nacme2_8=nacme

{ddr,2*dp
orbital,2140.2,2341.2,2342.2;
density,8000.2,8100.2,8200.2;
state, 4.1,3.1}
nacme2_9=nacme

{ddr,2*dp
orbital,2140.2,2341.2,2342.2;
density,8000.2,8100.2,8200.2;
state, 5.1,3.1}
nacme2_10=nacme

nacmeav7=(nacme7p+nacme7m)*0.5
nacmeav8=(nacme8p+nacme8m)*0.5
nacmeav9=(nacme9p+nacme9m)*0.5
nacmeav10=(nacme10p+nacme10m)*0.5

table,RHO,PHI,nacmeav7,nacmeav8,nacmeav9,nacmeav10,nacme2_7,nacme2_8,nacme2_9,nacme2_10
title,Non-adiabatic couplings for NO3
save,no3-nac-phi-q3-q4_2345.res

ENDDO
ENDDO
 ---
