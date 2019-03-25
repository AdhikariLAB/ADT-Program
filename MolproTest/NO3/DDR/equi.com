
***, Molpro template created from ADT program for "WhatEver"
memory,500,m
file,2,molpro.wfu,new;

basis=6-31G**;

symmetry,nosym


geometry=geom.xyz

{uhf;orbital,2200.2}
{multi;occ,19;closed,11; wf,31,1,1,0;state,5;orbital,2140.2;}
{mrci; occ,19;closed,11;core,11; wf,31,1,1,0;state,5;save,6000.2;noexc}
show,  energy
table, energy
save,equienr.res,new
{table,____; noprint,heading}

---

