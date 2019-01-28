import numpy as np


def jac2cart(rs,rc,gamma):
   """ Convert to cartesian coordinates in yz plane.
       Molecule will be on yz plane with m1 and m2 lying on y-axis.
       NOTE : coordinates will be in atomic units NOT angstroms """
   h1y = -rs/2.0
   h1z = 0.0
   h2y = rs/2.0
   h2z = 0.0
   h3y = rc*math.cos(gamma)
   h3z = rc*math.sin(gamma)
   lcart = [ (0,h1y,h1z), (0,h2y,h2z), (0,h3y,h3z) ]
   return lcart


def to_jacobi(rho,theta,phi):
    theta = np.deg2rad(theta)
    phi   = np.deg2rad(phi)
    mult  = np.sqrt(rho*rho/np.sqrt(3.0))

    R1 = np.sqrt((1.0+np.sin(theta)*np.sin(phi+4.0*np.pi/3.0)))*mult # H2-H3 distance
    R2 = np.sqrt((1.0+np.sin(theta)*np.sin(phi-4.0*np.pi/3.0)))*mult # H1-H3 distance
    R3 = np.sqrt((1.0+np.sin(theta)*np.sin(phi)))*mult               # H1-H2 distance


    if R3 > 1e-10:
        y3 = (R2*R2 + R3*R3 - R1*R1)/(2.0*R3)
    else:
        y3 = 0.0

    ss = R2*R2 - y3*y3

    if ss > -1e-10 and ss < 0.0:
        ss = 0.0

    z3 = np.sqrt(ss)

    rs = R3
    ycm = R3/2.0
    zcm = 0.0

    rc = np.sqrt((y3-ycm)**2 + (z3-zcm)**2)

    if rc < 1e-10:
        rc = 0.0

    if rs == 0.0 or rc == 0.0:
        gamma = 0.0
    else:
        val = (rc*rc+(rs*rs/4.0)-R1*R1)/(rc*rs)
        if val > 1.0 and val-1.0 < 1e-10:
            val = 1.0
        if val < -1.0 and abs(val)-1.0 < 1e-10:
            val = -1.0
        gamma=np.acos(val)

    return (rs, rc, gamma)