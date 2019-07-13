
***, Molpro template created from ADT program for analytical job
memory,10,m
file,2,molpro_init.wfu,new;

basis=6-31G*;

symmetry,x

geometry=geom.xyz
{hf}
{mcscf;occ,11,2;closed,3,0;wf,23,1,1,0;state,3;wf,23,2,1,0;state,2;;start,2140.2; orbital,2140.2;}


show, energy
table, energy


{table,____; noprint,heading}


---