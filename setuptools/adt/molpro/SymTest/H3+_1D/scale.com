
***, Molpro template created from ADT program for analytical job
memory,10,m
file,2,molpro_init.wfu,new;

basis=cc-pv5z;

symmetry,nosym

geometry=geom.xyz

{mcscf;occ,10;closed,0;;wf,2,1,0,1;state,3;;start,2140.2; orbital,2140.2;maxiter,40}
{mrci;occ,10;closed,0;;wf,2,1,0,1;state,3;;save,6000.2;dm,8000.2;}

show, energy
table, energy

save,scale.res,new
{table,____; noprint,heading}


---