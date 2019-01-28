***,C6H3F3 PES and non-adiabatic coupling terms
memory,100,m
file,2,c6h3f3_q9_2D_ana.wfu,new;
!gthresh,twoint=1.d-15,oneint=1.d-15,energy=1.d-10,gradient=1.d-6,prefac=1.d-15

basis=cc-pVDZ;

HCROSS=0.063508;
CMASS=12.011000;
HMASS=1.007940;
FMASS=18.998403;
CMINVTINV=0.001883651;

OMEGA9=1665*CMINVTINV;

X01=0.7040020239;
Y01=1.2193673894;
Z01=0.0000000000;
X02=1.3632629314;
Y02=0.0000000717;
Z02=0.0000000000;
X03=0.7040021010;
Y03=-1.2193673894;
Z03=0.0000000000;
X04=-0.6816314627;
Y04=-1.1806204519;
Z04=0.0000000000;
X05=-1.4080040150;
Y05=-0.0000000200;
Z05=0.0000000000;
X06=-0.6816315097;
Y06=1.1806203946;
Z06=0.0000000000;
X07=1.2435379779;
Y07=2.1538708923;
Z07=0.0000000000;
X08=1.2435378577;
Y08=-2.1538709952;
Z08=0.0000000000;
X09=-2.4870757626;
Y09=0.0000000470;
Z09=0.0000000000;
X010=2.7073025248;
Y010=0.0000000956;
Z010=0.0000000000;
X011=-1.3536513581;
Y011=-2.3445928579;
Z011=0.0000000000;
X012=-1.3536513086;
Y012=2.3445928238;
Z012=0.0000000000;

PROC CALCULATE_COORD

x1=SQRT(HCROSS/CMASS)*(-0.0834498*Q9X/SQRT(OMEGA9)-0.4481544*Q9Y/SQRT(OMEGA9));
y1=SQRT(HCROSS/CMASS)*( 0.2762947*Q9X/SQRT(OMEGA9)+0.1557110*Q9Y/SQRT(OMEGA9));
z1=SQRT(HCROSS/CMASS)*( 0.0000000*Q9X/SQRT(OMEGA9)-0.0000000*Q9Y/SQRT(OMEGA9));
x2=SQRT(HCROSS/CMASS)*(-0.0291331*Q9X/SQRT(OMEGA9)+0.2922013*Q9Y/SQRT(OMEGA9));
y2=SQRT(HCROSS/CMASS)*(-0.4900008*Q9X/SQRT(OMEGA9)-0.0488704*Q9Y/SQRT(OMEGA9));
z2=SQRT(HCROSS/CMASS)*(-0.0000000*Q9X/SQRT(OMEGA9)+0.0000000*Q9Y/SQRT(OMEGA9));
x3=SQRT(HCROSS/CMASS)*( 0.1703216*Q9X/SQRT(OMEGA9)-0.4228599*Q9Y/SQRT(OMEGA9));
y3=SQRT(HCROSS/CMASS)*( 0.3016045*Q9X/SQRT(OMEGA9)-0.0980774*Q9Y/SQRT(OMEGA9));
z3=SQRT(HCROSS/CMASS)*(-0.0000000*Q9X/SQRT(OMEGA9)-0.0000000*Q9Y/SQRT(OMEGA9));
x4=SQRT(HCROSS/CMASS)*( 0.0416696*Q9X/SQRT(OMEGA9)+0.4491374*Q9Y/SQRT(OMEGA9));
y4=SQRT(HCROSS/CMASS)*(-0.3331505*Q9X/SQRT(OMEGA9)-0.1197015*Q9Y/SQRT(OMEGA9));
z4=SQRT(HCROSS/CMASS)*(-0.0000000*Q9X/SQRT(OMEGA9)-0.0000000*Q9Y/SQRT(OMEGA9));
x5=SQRT(HCROSS/CMASS)*( 0.0215068*Q9X/SQRT(OMEGA9)-0.2156988*Q9Y/SQRT(OMEGA9));
y5=SQRT(HCROSS/CMASS)*( 0.5087999*Q9X/SQRT(OMEGA9)+0.0507429*Q9Y/SQRT(OMEGA9));
z5=SQRT(HCROSS/CMASS)*( 0.0000000*Q9X/SQRT(OMEGA9)-0.0000000*Q9Y/SQRT(OMEGA9));
x6=SQRT(HCROSS/CMASS)*(-0.1295525*Q9X/SQRT(OMEGA9)+0.4320494*Q9Y/SQRT(OMEGA9));
y6=SQRT(HCROSS/CMASS)*(-0.3502260*Q9X/SQRT(OMEGA9)+0.0515502*Q9Y/SQRT(OMEGA9));
z6=SQRT(HCROSS/CMASS)*( 0.0000000*Q9X/SQRT(OMEGA9)+0.0000000*Q9Y/SQRT(OMEGA9));

x7=SQRT(HCROSS/HMASS)*( 0.0863377*Q9X/SQRT(OMEGA9)+0.1017604*Q9Y/SQRT(OMEGA9));
y7=SQRT(HCROSS/HMASS)*( 0.0276167*Q9X/SQRT(OMEGA9)-0.0937339*Q9Y/SQRT(OMEGA9));
z7=SQRT(HCROSS/HMASS)*(-0.0000000*Q9X/SQRT(OMEGA9)+0.0000000*Q9Y/SQRT(OMEGA9));
x8=SQRT(HCROSS/HMASS)*(-0.1047342*Q9X/SQRT(OMEGA9)+0.0827080*Q9Y/SQRT(OMEGA9));
y8=SQRT(HCROSS/HMASS)*( 0.0085591*Q9X/SQRT(OMEGA9)+0.0973467*Q9Y/SQRT(OMEGA9));
z8=SQRT(HCROSS/HMASS)*( 0.0000000*Q9X/SQRT(OMEGA9)+0.0000000*Q9Y/SQRT(OMEGA9));
x9=SQRT(HCROSS/HMASS)*( 0.0073059*Q9X/SQRT(OMEGA9)-0.0732420*Q9Y/SQRT(OMEGA9));
y9=SQRT(HCROSS/HMASS)*(-0.1474979*Q9X/SQRT(OMEGA9)-0.0147085*Q9Y/SQRT(OMEGA9));
z9=SQRT(HCROSS/HMASS)*( 0.0000000*Q9X/SQRT(OMEGA9)-0.0000000*Q9Y/SQRT(OMEGA9));

x10=SQRT(HCROSS/FMASS)*( 0.0039163*Q9X/SQRT(OMEGA9)-0.0393580*Q9Y/SQRT(OMEGA9));
y10=SQRT(HCROSS/FMASS)*( 0.0236702*Q9X/SQRT(OMEGA9)+0.0023617*Q9Y/SQRT(OMEGA9));
z10=SQRT(HCROSS/FMASS)*(-0.0000000*Q9X/SQRT(OMEGA9)+0.0000000*Q9Y/SQRT(OMEGA9));
x11=SQRT(HCROSS/FMASS)*( 0.0095630*Q9X/SQRT(OMEGA9)-0.0269101*Q9Y/SQRT(OMEGA9));
y11=SQRT(HCROSS/FMASS)*( 0.0361245*Q9X/SQRT(OMEGA9)-0.0032703*Q9Y/SQRT(OMEGA9));
z11=SQRT(HCROSS/FMASS)*( 0.0000000*Q9X/SQRT(OMEGA9)+0.0000000*Q9Y/SQRT(OMEGA9));
x12=SQRT(HCROSS/FMASS)*(-0.0040570*Q9X/SQRT(OMEGA9)-0.0282680*Q9Y/SQRT(OMEGA9));
y12=SQRT(HCROSS/FMASS)*( 0.0347659*Q9X/SQRT(OMEGA9)+0.0103382*Q9Y/SQRT(OMEGA9));
z12=SQRT(HCROSS/FMASS)*(-0.0000000*Q9X/SQRT(OMEGA9)+0.0000000*Q9Y/SQRT(OMEGA9));

ENDPROC

!------------------------------------------------------------
!EQUILIBRIUM POINT
!------------------------------------------------------------

Q9X=0.0
Q9Y=0.0

CALCULATE_COORD

symmetry,nosym

geomtyp=xyz
geometry= {
12
This is the geometry input for C6H3F3 in XYZ format
c,x1+x01,y1+y01,z1+z01
c,x2+x02,y2+y02,z2+z02
c,x3+x03,y3+y03,z3+z03
c,x4+x04,y4+y04,z4+z04
c,x5+x05,y5+y05,z5+z05
c,x6+x06,y6+y06,z6+z06
h,x7+x07,y7+y07,z7+z07
h,x8+x08,y8+y08,z8+z08
h,x9+x09,y9+y09,z9+z09
f,x10+x010,y10+y010,z10+z010
f,x11+x011,y11+y011,z11+z011
f,x12+x012,y12+y012,z12+z012
}

{uhf;wf,65,1,1,1;accu,5;orbital,2200.2}
{multi;occ,36;closed,27;wf,65,1,1,1;state,6;
       start,2200.2;orbital,2340.2;maxiter,40;}

!show, energy
!table, energy
!save,c6h3f3-enrq9-mcscf.res
!title 'results: EQUILIBRIUM'
!{table,____; noprint,heading}

{mrci;occ,33;core,27;closed,27;wf,65,1,1,1;state,6;}

!show, energy
!table, energy
!save,c6h3f3-enrq9-mrci.res
!title 'results: EQUILIBRIUM'
!{table,____; noprint,heading}

!------------------------------------------------------------
!GRIDS
!------------------------------------------------------------

RHOMIN=0.0
RHOMAX=1.0
DRHO=0.02
PHIMIN=0.0
PHIMAX=360.0
DPHI=2.0

DO I=1,2
PHI=PHIMIN+I*DPHI
DO J=1,50
RHO=RHOMIN+J*DRHO

Q9X=RHO*COS(PHI);
Q9Y=RHO*SIN(PHI);

CALCULATE_COORD

symmetry,nosym

geomtyp=xyz
geometry= {
12
This is the geometry input for C6H3F3 in XYZ format
c,x1+x01,y1+y01,z1+z01
c,x2+x02,y2+y02,z2+z02
c,x3+x03,y3+y03,z3+z03
c,x4+x04,y4+y04,z4+z04
c,x5+x05,y5+y05,z5+z05
c,x6+x06,y6+y06,z6+z06
h,x7+x07,y7+y07,z7+z07
h,x8+x08,y8+y08,z8+z08
h,x9+x09,y9+y09,z9+z09
f,x10+x010,y10+y010,z10+z010
f,x11+x011,y11+y011,z11+z011
f,x12+x012,y12+y012,z12+z012
}

basis=cc-pVDZ;

if(j.eq.1)then
{multi;occ,36;closed,27;wf,65,1,1,1;state,6;
       start,2340.2;orbital,2140.2+j;maxiter,40;}
else
{multi;occ,36;closed,27;wf,65,1,1,1;state,6;
       start,2140.2+(j-1);orbital,2140.2+j;maxiter,40;}
endif

show, energy
table, energy
save,c6h3f3-enrq9-mcscf.res
title 'results: Geom --- RHO= $RHO PHI= $PHI
{table,____; noprint,heading}

{mrci;occ,33;core,27;closed,27;wf,65,1,1,1;state,6;MAXITER,99,99;
             save,6000.2;dm,8000.2}

show, energy
table, energy
save,c6h3f3-enrq9-mrci.res
title 'results: Geom --- RHO= $RHO PHI= $PHI
{table,____; noprint,heading}

basis=6-31G**

{multi;occ,36;closed,27;wf,65,1,1,1;state,6;start,2140.2+j;orbital,2350.2;maxiter,40;
cpmcscf,nacm,1.1,2.1,accu=1.d-10,record=5101.1;
cpmcscf,nacm,1.1,3.1,accu=1.d-10,record=5102.1;
cpmcscf,nacm,1.1,4.1,accu=1.d-10,record=5103.1;
cpmcscf,nacm,1.1,5.1,accu=1.d-10,record=5104.1;
cpmcscf,nacm,1.1,6.1,accu=1.d-10,record=5105.1;
}

{multi;occ,36;closed,27;wf,65,1,1,1;state,6;start,2140.2+j;orbital,2350.2;maxiter,40;
cpmcscf,nacm,2.1,3.1,accu=1.d-10,record=5106.1;
cpmcscf,nacm,2.1,4.1,accu=1.d-10,record=5107.1;
cpmcscf,nacm,2.1,5.1,accu=1.d-10,record=5108.1;
cpmcscf,nacm,2.1,6.1,accu=1.d-10,record=5109.1;
cpmcscf,nacm,3.1,4.1,accu=1.d-10,record=5110.1;
}

{multi;occ,36;closed,27;wf,65,1,1,1;state,6;start,2140.2+j;orbital,2350.2;maxiter,40;
cpmcscf,nacm,3.1,5.1,accu=1.d-10,record=5111.1;
cpmcscf,nacm,3.1,6.1,accu=1.d-10,record=5112.1;
cpmcscf,nacm,4.1,5.1,accu=1.d-10,record=5113.1;
cpmcscf,nacm,4.1,6.1,accu=1.d-10,record=5114.1;
cpmcscf,nacm,5.1,6.1,accu=1.d-10,record=5115.1;
}

force;nacm,5101.1;varsav
table,gradx,grady,gradz;
save,c6h3f3-analytic-nac.res;
title 'results: Tau_12 ---Geom --- RHO= $RHO PHI= $PHI
{table,____;  noprint,heading}

force;nacm,5102.1;varsav
table,gradx,grady,gradz;
save,c6h3f3-analytic-nac.res;
title 'results: Tau_13 ---Geom --- RHO= $RHO PHI= $PHI
{table,____; noprint,heading}

force;nacm,5103.1;varsav
table,gradx,grady,gradz;
save,c6h3f3-analytic-nac.res;
title 'results: Tau_14 ---Geom --- RHO= $RHO PHI= $PHI
{table,____; noprint,heading}

force;nacm,5104.1;varsav
table,gradx,grady,gradz;
save,c6h3f3-analytic-nac.res;
title 'results: Tau_15 ---Geom --- RHO= $RHO PHI= $PHI
{table,____; noprint,heading}

force;nacm,5105.1;varsav
table,gradx,grady,gradz;
save,c6h3f3-analytic-nac.res;
title 'results: Tau_16 ---Geom --- RHO= $RHO PHI= $PHI
{table,____; noprint,heading}

force;nacm,5106.1;varsav
table,gradx,grady,gradz;
save,c6h3f3-analytic-nac.res;
title 'results: Tau_23 ---Geom --- RHO= $RHO PHI= $PHI
{table,____; noprint,heading}

force;nacm,5107.1;varsav
table,gradx,grady,gradz;
save,c6h3f3-analytic-nac.res;
title 'results: Tau_24 ---Geom --- RHO= $RHO PHI= $PHI
{table,____; noprint,heading}

force;nacm,5108.1;varsav
table,gradx,grady,gradz;
save,c6h3f3-analytic-nac.res;
title 'results: Tau_25 ---Geom --- RHO= $RHO PHI= $PHI
{table,____; noprint,heading}

force;nacm,5109.1;varsav
table,gradx,grady,gradz;
save,c6h3f3-analytic-nac.res;
title 'results: Tau_26 ---Geom --- RHO= $RHO PHI= $PHI
{table,____; noprint,heading}

force;nacm,5110.1;varsav
table,gradx,grady,gradz;
save,c6h3f3-analytic-nac.res;
title 'results: Tau_34 ---Geom --- RHO= $RHO PHI= $PHI
{table,____; noprint,heading}

force;nacm,5111.1;varsav
table,gradx,grady,gradz;
save,c6h3f3-analytic-nac.res;
title 'results: Tau_35 ---Geom --- RHO= $RHO PHI= $PHI
{table,____; noprint,heading}

force;nacm,5112.1;varsav
table,gradx,grady,gradz;
save,c6h3f3-analytic-nac.res;
title 'results: Tau_36 ---Geom --- RHO= $RHO PHI= $PHI
{table,____; noprint,heading}

force;nacm,5113.1;varsav
table,gradx,grady,gradz;
save,c6h3f3-analytic-nac.res;
title 'results: Tau_45 ---Geom --- RHO= $RHO PHI= $PHI
{table,____; noprint,heading}

force;nacm,5114.1;varsav
table,gradx,grady,gradz;
save,c6h3f3-analytic-nac.res;
title 'results: Tau_46 ---Geom --- RHO= $RHO PHI= $PHI
{table,____; noprint,heading}

force;nacm,5115.1;varsav
table,gradx,grady,gradz;
save,c6h3f3-analytic-nac.res;
title 'results: Tau_56 ---Geom --- RHO= $RHO PHI= $PHI
{table,____; noprint,heading}

ENDDO
ENDDO
 ---
