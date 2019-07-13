
***, Molpro template created from ADT program for analytical job
memory,10,m
file,2,molpro_init.wfu,new;

basis=cc-pv5z;

symmetry,x

geometry=geom.xyz

{mcscf;occ,10,1;wf,2,1,0,1;state,3;;start,2140.2; orbital,2140.2;maxiter,40;ACCURACY,GRADIENT=1.d-4,ENERGY=1.d-6,STEP=1.d-2}
        {mrci;occ,10,1;wf,2,1,0,1;state,3;;save,6000.2;dm,8000.2;}


show, energy
table, energy
save,enr.res,new

{table,____; noprint,heading}



basis={spd,h,cc-pv5z}

{mcscf;occ,10,1; wf,2,1,0,1;state,3;; start,2140.2;maxiter,40;ACCURACY,GRADIENT=1.d-4,ENERGY=1.d-6,STEP=1.d-2;
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


---
