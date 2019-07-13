
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
save,enr.res,new

{table,____; noprint,heading}



basis=cc-pv5z


!for +d1
symmetry,nosym
geometry=geom2.xyz
{multi;occ,10;closed,0; ;wf,2,1,0,1;state,3;;start,2140.2;orbital,2241.2;maxiter,40}
{mrci; occ,10;closed,0;; wf,2,1,0,1;state,3;;save,6001.2;}
{ci;trans,6000.2,6001.2;dm,8001.2;}


!for -d1
symmetry,nosym
geometry=geom3.xyz
{multi;occ,10;closed,0;; wf,2,1,0,1;state,3;;start,2140.2;orbital,2242.2;maxiter,40}
{mrci; occ,10;closed,0;; wf,2,1,0,1;state,3;;save,6002.2;}
{ci;trans,6000.2,6002.2;dm,8002.2;}


!N1D:
!for +d2
symmetry,nosym
geometry=geom4.xyz
{multi;occ,10;closed,0;; wf,2,1,0,1;state,3;;start,2140.2;orbital,2243.2;maxiter,40}
{mrci; occ,10;closed,0;; wf,2,1,0,1;state,3;;save,6003.2;}
{ci;trans,6000.2,6003.2;dm,8003.2;}


!N1D:
!for -d2
symmetry,nosym
geometry=geom5.xyz
{multi;occ,10;closed,0;; wf,2,1,0,1;state,3;;start,2140.2;orbital,2244.2;maxiter,40}
{mrci; occ,10;closed,0;; wf,2,1,0,1;state,3;;save,6004.2;}
{ci;trans,6000.2,6004.2;dm,8004.2;}



!for taur     
{ddr, 2*0.01
orbital,2140.2,2141.2,2142.2;
density,8000.2,8001.2,8002.2;
state, 2.1,1.1
}
nacmr = nacme


!N1D:for taup
{ddr, 2*0.03
orbital,2140.2,2143.2,2144.2;
density,8000.2,8003.2,8004.2;
state, 2.1,1.1
}
nacmp = nacme

table, nacmr,nacmp
save,ddrnact_1.1_2.1.res,new;


!for taur     
{ddr, 2*0.01
orbital,2140.2,2141.2,2142.2;
density,8000.2,8001.2,8002.2;
state, 3.1,1.1
}
nacmr = nacme


!N1D:for taup
{ddr, 2*0.03
orbital,2140.2,2143.2,2144.2;
density,8000.2,8003.2,8004.2;
state, 3.1,1.1
}
nacmp = nacme

table, nacmr,nacmp
save,ddrnact_1.1_3.1.res,new;


!for taur     
{ddr, 2*0.01
orbital,2140.2,2141.2,2142.2;
density,8000.2,8001.2,8002.2;
state, 3.1,2.1
}
nacmr = nacme


!N1D:for taup
{ddr, 2*0.03
orbital,2140.2,2143.2,2144.2;
density,8000.2,8003.2,8004.2;
state, 3.1,2.1
}
nacmp = nacme

table, nacmr,nacmp
save,ddrnact_2.1_3.1.res,new;


---
