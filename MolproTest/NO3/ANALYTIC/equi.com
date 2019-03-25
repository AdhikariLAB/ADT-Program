
***, Molpro template created from ADT program for "WhatEver"
file,2,molpro.wfu,new;

basis=6-31G**;

symmetry,nosym


geometry=geom.xyz

{uhf;orbital,2200.2}
{multi;occ,19;closed,11; wf,31,1,1,0;state,5;orbital,2140.2;}

show, energy
table, energy
save,equienr.res,new
{table,____; noprint,heading}

---

