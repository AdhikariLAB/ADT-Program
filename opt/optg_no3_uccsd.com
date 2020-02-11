!Geometry optimization of NO3 using basis set cc-pVTZ 
memory,1000,m



BASIS=cc-pVTZ;

symmetry,nosym

geometry={
N
O1,N,rno1
O2,N,rno2,O1,hno21
O3,N,rno3,O1,hno31,O2,180
}

rno1=1.24 angstrom
rno2=1.24 angstrom
rno3=1.24 angstrom
hno21=120 degree
hno31=120 degree

uhf;
uccsd(t):ORBITAL,2200.1;
{optg,savexyz=optg.xyz}
{frequencies;print,hessian}

---
