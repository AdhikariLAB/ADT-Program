
***, Molpro template created from ADT program for "WhatEver"
memory,100,m
file,2,molpro.wfu;

basis=6-31G**;

symmetry,nosym


geometry=geom.xyz

{multi;occ,19;closed,11; wf,31,1,1,0;state,5;start,2140.2; orbital,2140.2;}

show, energy
table, energy
save,enr.res,new
{table,____; noprint,heading}


basis=6-31G**

{multi;occ,19;closed,11; wf,31,1,1,0;state,5; start,2140.2;
cpmcscf,nacm,1.1,2.1,record=5112.1;
cpmcscf,nacm,1.1,3.1,record=5113.1;
cpmcscf,nacm,1.1,4.1,record=5114.1;
cpmcscf,nacm,1.1,5.1,record=5115.1;
cpmcscf,nacm,2.1,3.1,record=5123.1;
}

force;nacm,5112.1;varsav
table,gradx,grady,gradz;
save,ananac12.res,new;
{table,____;  noprint,heading}


force;nacm,5113.1;varsav
table,gradx,grady,gradz;
save,ananac13.res,new;
{table,____;  noprint,heading}


force;nacm,5114.1;varsav
table,gradx,grady,gradz;
save,ananac14.res,new;
{table,____;  noprint,heading}


force;nacm,5115.1;varsav
table,gradx,grady,gradz;
save,ananac15.res,new;
{table,____;  noprint,heading}


force;nacm,5123.1;varsav
table,gradx,grady,gradz;
save,ananac23.res,new;
{table,____;  noprint,heading}


{multi;occ,19;closed,11; wf,31,1,1,0;state,5; start,2140.2;
cpmcscf,nacm,2.1,4.1,record=5124.1;
cpmcscf,nacm,2.1,5.1,record=5125.1;
cpmcscf,nacm,3.1,4.1,record=5134.1;
cpmcscf,nacm,3.1,5.1,record=5135.1;
cpmcscf,nacm,4.1,5.1,record=5145.1;
}

force;nacm,5124.1;varsav
table,gradx,grady,gradz;
save,ananac24.res,new;
{table,____;  noprint,heading}


force;nacm,5125.1;varsav
table,gradx,grady,gradz;
save,ananac25.res,new;
{table,____;  noprint,heading}


force;nacm,5134.1;varsav
table,gradx,grady,gradz;
save,ananac34.res,new;
{table,____;  noprint,heading}


force;nacm,5135.1;varsav
table,gradx,grady,gradz;
save,ananac35.res,new;
{table,____;  noprint,heading}


force;nacm,5145.1;varsav
table,gradx,grady,gradz;
save,ananac45.res,new;
{table,____;  noprint,heading}


---
