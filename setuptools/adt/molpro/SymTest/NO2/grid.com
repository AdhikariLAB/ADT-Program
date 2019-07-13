
***, Molpro template created from ADT program for analytical job
memory,10,m
file,2,molpro_init.wfu,new;

basis=6-31G*;

symmetry,x

geometry=geom.xyz

{mcscf;occ,11,2;closed,3,0;wf,23,1,1,0;state,3;wf,23,2,1,0;state,2;;start,2140.2; orbital,2140.2;}


show, energy
table, energy
save,enr.res,new

{table,____; noprint,heading}



basis=6-31G*

{mcscf;occ,11,2;closed,3,0; wf,23,1,1,0;state,3;wf,23,2,1,0;state,2;; start,2140.2;;
cpmcscf,nacm,1.1,2.1,record=5101.1,;
cpmcscf,nacm,1.1,3.1,record=5102.1,;
cpmcscf,nacm,2.1,3.1,record=5103.1,;
}


force;nacm,5101.1;varsav
table,gradx,grady,gradz;
save,ananac_1.1_2.1.res,new;
{table,____;  noprint,heading}


force;nacm,5102.1;varsav
table,gradx,grady,gradz;
save,ananac_1.1_3.1.res,new;
{table,____;  noprint,heading}


force;nacm,5103.1;varsav
table,gradx,grady,gradz;
save,ananac_2.1_3.1.res,new;
{table,____;  noprint,heading}


{mcscf;occ,11,2;closed,3,0; wf,23,1,1,0;state,3;wf,23,2,1,0;state,2;; start,2140.2;;
cpmcscf,nacm,1.2,2.2,record=5101.1,;
}


force;nacm,5101.1;varsav
table,gradx,grady,gradz;
save,ananac_1.2_2.2.res,new;
{table,____;  noprint,heading}


---
