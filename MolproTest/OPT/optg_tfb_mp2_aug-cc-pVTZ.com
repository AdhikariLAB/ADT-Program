!Geometry optimization of C6H3F3 using basis set aug-cc-pvTZ
memory,100,m

file,2,hess_c6h3f3_mp2.wfu,new;

GPRINT,DISTANCE,ANGLES,ORBITAL,ORBEN,CIVECTOR,REF;


symmetry, nosym

geometry={
N
O1,N,rno1
O2,N,rno2,O1,alph
}

rno1 = 1.2126 ANGSTROM
rno2 = 1.2126 ANGSTROM
alph = 133.7 DEGREE


hf

optg;coord,3n

{frequencies
print,hessian}

---
