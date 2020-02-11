memory,1000,m

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

{ks,b3lyp}
{optg,savexyz=optg.xyz}
{frequencies;print,hessian}

---
