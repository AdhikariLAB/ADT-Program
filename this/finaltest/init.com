
***, Molpro template created from ADT program
memory,10,m
file,2,molpro_init.wfu;

basis=6-31G**;

symmetry,nosym


geometry=geom.xyz
{uhf}
{multi;occ,19;closed,11; wf,31,1,1,0;state,5;start,2140.2; orbital,2140.2;}

show, energy
table, energy
save,enr.res,new
{table,____; noprint,heading}


basis=6-31G**

{mcscf;occ,19;closed,11; wf,31,1,1,0;state,5; start,2140.2;
cpmcscf,nacm,1.1,2.1,record=5101.1;
cpmcscf,nacm,1.1,3.1,record=5102.1;
cpmcscf,nacm,1.1,4.1,record=5103.1;
cpmcscf,nacm,1.1,5.1,record=5104.1;
cpmcscf,nacm,2.1,3.1,record=5105.1;
}

force;nacm,5101.1;varsav
table,gradx,grady,gradz;
save,ananac12.res,new;
{table,____;  noprint,heading}


force;nacm,5102.1;varsav
table,gradx,grady,gradz;
save,ananac13.res,new;
{table,____;  noprint,heading}


force;nacm,5103.1;varsav
table,gradx,grady,gradz;
save,ananac14.res,new;
{table,____;  noprint,heading}


force;nacm,5104.1;varsav
table,gradx,grady,gradz;
save,ananac15.res,new;
{table,____;  noprint,heading}


force;nacm,5105.1;varsav
table,gradx,grady,gradz;
save,ananac23.res,new;
{table,____;  noprint,heading}


{mcscf;occ,19;closed,11; wf,31,1,1,0;state,5; start,2140.2;
cpmcscf,nacm,2.1,4.1,record=5101.1;
cpmcscf,nacm,2.1,5.1,record=5102.1;
cpmcscf,nacm,3.1,4.1,record=5103.1;
cpmcscf,nacm,3.1,5.1,record=5104.1;
cpmcscf,nacm,4.1,5.1,record=5105.1;
}

force;nacm,5101.1;varsav
table,gradx,grady,gradz;
save,ananac24.res,new;
{table,____;  noprint,heading}


force;nacm,5102.1;varsav
table,gradx,grady,gradz;
save,ananac25.res,new;
{table,____;  noprint,heading}


force;nacm,5103.1;varsav
table,gradx,grady,gradz;
save,ananac34.res,new;
{table,____;  noprint,heading}


force;nacm,5104.1;varsav
table,gradx,grady,gradz;
save,ananac35.res,new;
{table,____;  noprint,heading}


force;nacm,5105.1;varsav
table,gradx,grady,gradz;
save,ananac45.res,new;
{table,____;  noprint,heading}


---
