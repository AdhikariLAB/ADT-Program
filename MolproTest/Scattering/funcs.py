import numpy as np


def to_jacobi(rho, theta, phi, m1, m2, m3):
    """ returns jacobi coordinates """

    M = m1 + m2 + m3  #use global variables of m1,m2,m3
    mu = np.sqrt(m1*m2*m3/M)
    d1 = np.sqrt(m1*(m2+m3)/(mu*M))
    d2 = np.sqrt(m2*(m3+m1)/(mu*M))
    d3 = np.sqrt(m3*(m1+m2)/(mu*M))
    eps3 = 2*np.atan(m2/mu)
    eps2 = 2*np.atan(m3/mu)

    R1 = (1.0/np.sqrt(2.0))*rho*d3*np.sqrt(1.0+np.sin(theta)*np.cos(phi+eps3))
    R2 = (1.0/np.sqrt(2.0))*rho*d1*np.sqrt(1.0+np.sin(theta)*np.cos(phi))      
    R3 = (1.0/np.sqrt(2.0))*rho*d2*np.sqrt(1.0+np.sin(theta)*np.cos(phi-eps2)) 

    if R1 < 1e-10:
        R1 = 0.0
    if R2 < 1e-10:
        R2 = 0.0
    if R3 < 1e-10:
        R3 = 0.0

    area = AreaTriangle(R1,R2,R3)
    x = R2*R2 + R3*R3 - R1*R1
    y = 4.0*area
    Ang123 = np.atan2(y,x)
    x2 = (0.0,0.0)
    x3 = (R2,0.0)
    x1 = (R3*np.cos(Ang123),R3*np.sin(Ang123))
    # these are non-mass scaled jacobi coords
    # r : (x3-x2)
    # R : x1 - com(x3,x2)
    # gamma : angle between r and R
    r = (R2,0.0)
    R = (R3*np.cos(Ang123) - m3*R2/(m2+m3) , R3*np.sin(Ang123))
    rs = np.sqrt(dotp(r,r))
    rc = np.sqrt(dotp(R,R))
    if rc < 1e-10:
        rc = 0.0

    rtil = (R2*m2/(m2+m3),0.0)
    drtil = np.sqrt(dotp(rtil,rtil))
    Areasmall = AreaTriangle(drtil,R1,rc)
    y = 4.0 * Areasmall
    x = drtil*drtil + rc*rc - R1*R1
    if np.fabs(x) < 1e-10:
        x = 0.0
    gamma = np.atan2(y,x)

    # we assert that gamma is always less than 90 degree always - our grid does not produce this for H2F
    # assert (gamma <= np.pi/2)

    return (rs, rc, gamma)


def jacobixyz(rs,rc,gamm):  # rs is jacobi 'r', rc is jacobi 'R' gamm is jacobi 'gamma'
# xyz coordinates generated from jacobi coordinates
# use global variables m1,m2,m3 to assign the masses of the atoms
p1 = m1,(rc*np.cos(gamm),rc*np.sin(gamm))
p2 = m2,( -rs/2.0 , 0.0 )
p3 = m3,( rs/2.0 , 0.0 )
return p1,p2,p3
