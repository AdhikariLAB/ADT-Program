
***, Molpro template created from ADT program for "WhatEver"
memory,10,m
file,2,molpro.wfu;

basis=6-311++G**;

symmetry,nosym

geomtyp=xyz
geometry=geom1.xyz

!replacebyuhf
{multi;occ,11;closed,1; wf,11,1,1,0;state,3;orbital,2140.2;}
{mrci; occ,11;closed,1; wf,11,1,1,0;state,3;save,6000.2;noexc}

show, energy
table, energy
save,enr.res,new
{table,____; noprint,heading}

!for +dr
symmetry,nosym
geometry=geom2.xyz
{multi;occ,11;closed,1 ;wf,11,1,1,0;state,3;start,2140.2;orbital,2241.2;}
{mrci; occ,11;closed,1; wf,11,1,1,0;state,3;save,6001.2;noexc}
{ci;trans,6000.2,6001.2;dm,8001.2}


!for -dr
symmetry,nosym
geometry=geom3.xyz
{multi;occ,11;closed,1; wf,11,1,1,0;state,3;start,2140.2;orbital,2242.2;}
{mrci; occ,11;closed,1; wf,11,1,1,0;state,3;save,6002.2;noexc}
{ci;trans,6000.2,6002.2;dm,8002.2}



!for +dp
symmetry,nosym
geometry=geom4.xyz
{multi;occ,11;closed,1; wf,11,1,1,0;state,3;start,2140.2;orbital,2243.2;}
{mrci; occ,11;closed,1; wf,11,1,1,0;state,3;save,6003.2;noexc}
{ci;trans,6000.2,6003.2;dm,8003.2}


!for -dp
symmetry,nosym
geometry=geom5.xyz
{multi;occ,11;closed,1; wf,11,1,1,0;state,3;start,2140.2;orbital,2244.2;}
{mrci; occ,11;closed,1; wf,11,1,1,0;state,3;save,6004.2;noexc}
{ci;trans,6000.2,6004.2;dm,8004.2}



!for taur     
{ddr,0.01,2140.2,2241.2,8001.2;state, 2.1,1.1}
nacmepv=nacme

{ddr,-0.01,2140.2,2242.2,8002.2;state, 2.1,1.1}
nacmemv=nacme

nacmr = 0.5*(nacmepv+ nacmemv)

!for taup
{ddr,0.02,2140.2,2243.2,8003.2;state, 2.1,1.1}
nacmepv=nacme

{ddr,-0.02,2140.2,2244.2,8004.2;state, 2.1,1.1}
nacmemv=nacme
nacmp = 0.5*(nacmepv+ nacmemv)


table, nacmr,nacmp
save,ddrnact12.res,new;



!for taur     
{ddr,0.01,2140.2,2241.2,8001.2;state, 3.1,1.1}
nacmepv=nacme

{ddr,-0.01,2140.2,2242.2,8002.2;state, 3.1,1.1}
nacmemv=nacme

nacmr = 0.5*(nacmepv+ nacmemv)

!for taup
{ddr,0.02,2140.2,2243.2,8003.2;state, 3.1,1.1}
nacmepv=nacme

{ddr,-0.02,2140.2,2244.2,8004.2;state, 3.1,1.1}
nacmemv=nacme
nacmp = 0.5*(nacmepv+ nacmemv)


table, nacmr,nacmp
save,ddrnact13.res,new;



!for taur     
{ddr,0.01,2140.2,2241.2,8001.2;state, 3.1,2.1}
nacmepv=nacme

{ddr,-0.01,2140.2,2242.2,8002.2;state, 3.1,2.1}
nacmemv=nacme

nacmr = 0.5*(nacmepv+ nacmemv)

!for taup
{ddr,0.02,2140.2,2243.2,8003.2;state, 3.1,2.1}
nacmepv=nacme

{ddr,-0.02,2140.2,2244.2,8004.2;state, 3.1,2.1}
nacmemv=nacme
nacmp = 0.5*(nacmepv+ nacmemv)


table, nacmr,nacmp
save,ddrnact23.res,new;



---
