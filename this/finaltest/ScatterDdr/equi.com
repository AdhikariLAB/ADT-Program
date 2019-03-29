
***, Molpro template created from ADT program
memory,10,m
file,2,molpro_equi.wfu,new;

basis=6-311++G**;

symmetry,nosym


geometry=geom.xyz

{uhf;orbital,2200.2}
{multi;occ,11;closed,1; wf,11,1,1,0;state,3;orbital,2140.2;}
{mrci; occ,11;closed,1; wf,11,1,1,0;state,3;save,6000.2;noexc}
show,  energy
table, energy
save,equienr.res,new
{table,____; noprint,heading}

---

