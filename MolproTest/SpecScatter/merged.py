import os
import sys
import shutil
import textwrap
import subprocess
import numpy as np
import ConfigParser



class Base():
    angtobohr = 1.8897259886
    hbar = 0.063508
    cInvToTauInv = 0.001883651

    def sin(self, x):
        """ A sin function that directly takes degree as input unlike numpy"""
        return np.sin(np.deg2rad(x))

    def cos(self, x):
        """ A cos function that directly takes degree as input unlike numpy"""
        return np.cos(np.deg2rad(x))

    def createAnaTemplate(self):
        molproEquiTemplate=textwrap.dedent('''
            ***, Molpro template created from ADT program
            memory,{memory};
            file,2,molpro_equi.wfu,new;

            basis={basis};

            symmetry,nosym

            geometry=geom.xyz

            {{uhf;orbital,2200.2}}
            {{{method};{cas}; wf,{wf};state,{state};orbital,2140.2;}}

            show, energy
            table, energy
            save,equienr.res,new
            {{table,____; noprint,heading}}

            ---

            '''.format(memory = self.memory,
                        basis = self.eInfo['basis'],
                        method=self.eInfo['method'],
                        state=self.eInfo['state'],
                        wf =  self.eInfo['wf'],
                        cas = self.eInfo['cas']))

        with open('equi.com' ,'w') as f:
            f.write(molproEquiTemplate)



        molproGridTemplate=textwrap.dedent('''
            ***, Molpro template created from ADT program
            memory,{memory}
            file,2,molpro.wfu;

            basis={basis};

            symmetry,nosym


            geometry=geom.xyz

            {{{method};{cas}; wf,{wf};state,{state};start,2140.2; orbital,2140.2;}}

            show, energy
            table, energy
            save,enr.res,new
            {{table,____; noprint,heading}}
            '''.format(memory = self.memory,
                        basis = self.eInfo['basis'],
                        method=self.eInfo['method'],
                        state=self.eInfo['state'],
                        wf =  self.eInfo['wf'],
                        cas = self.eInfo['cas']))


        nactTemp= "\n\nbasis={}\n".format(self.nInfo['basis'])

        for ind in range(0,len(self.nactPairs),5):
            nactTemp += "\n{{{method};{cas}; wf,{wf};state,{state}; start,2140.2;\n".format(
                                method=self.nInfo['method'],
                                state=self.eInfo['state'],
                                wf =  self.eInfo['wf'],
                                cas = self.eInfo['cas'])

            pairChunk = self.nactPairs[ind:ind+5]
            forceTemp = ''
            for count,pair in enumerate(pairChunk, start=1):
                f,s = pair
                nactTemp += "{nmethod},nacm,{f}.1,{s}.1,record=51{n:02}.1;\n".format(
                            nmethod = self.nInfo['method'],f=f, s=s, n=count)

                forceTemp +=textwrap.dedent("""
                force;nacm,51{n:02}.1;varsav
                table,gradx,grady,gradz;
                save,ananac{f}{s}.res,new;
                {{table,____;  noprint,heading}}

                """.format(f=f, s=s, n=count))
            nactTemp += "}\n" + forceTemp




        molproGridTemplate += nactTemp + '\n---\n'

        with open('grid.com','w') as f:
            f.write(molproGridTemplate)

    #which basis for which nact ????
    def createDdrTemplate(self):
        molproEquiTemplate=textwrap.dedent('''
            ***, Molpro template created from ADT program
            memory,{memory}
            file,2,molpro_equi.wfu,new;

            basis={basis};

            symmetry,nosym


            geometry=geom.xyz

            {{uhf;orbital,2200.2}}
            {{multi;{cas}; wf,{wf};state,{state};orbital,2140.2;}}
            {{mrci; {cas}; wf,{wf};state,{state};save,6000.2;noexc}}
            show,  energy
            table, energy
            save,equienr.res,new
            {{table,____; noprint,heading}}

            ---

            '''.format( memory = self.memory,
                        basis = self.eInfo['basis'],
                        state=self.eInfo['state'],
                        wf =  self.eInfo['wf'],
                        cas = self.eInfo['cas']))

        with open('equi.com' ,'w') as f:
            f.write(molproEquiTemplate)



        molproGridTemplate=textwrap.dedent('''
            ***, Molpro template created from ADT program for "WhatEver"
            memory,{memory}
            file,2,molpro.wfu;

            basis={basis};

            symmetry,nosym

            geomtyp=xyz
            geometry=geom1.xyz

            {{multi;{cas}; wf,{wf};state,{state};orbital,2140.2;}}
            {{mrci; {cas}; wf,{wf};state,{state};save,6000.2;noexc}}

            show, energy
            table, energy
            save,enr.res,new
            {{table,____; noprint,heading}}

            !for +dr
            symmetry,nosym
            geometry=geom2.xyz
            {{multi;{cas} ;wf,{wf};state,{state};start,2140.2;orbital,2241.2;}}
            {{mrci; {cas}; wf,{wf};state,{state};save,6001.2;noexc}}
            {{ci;trans,6000.2,6001.2;dm,8001.2}}


            !for -dr
            symmetry,nosym
            geometry=geom3.xyz
            {{multi;{cas}; wf,{wf};state,{state};start,2140.2;orbital,2242.2;}}
            {{mrci; {cas}; wf,{wf};state,{state};save,6002.2;noexc}}
            {{ci;trans,6000.2,6002.2;dm,8002.2}}



            !for +dp
            symmetry,nosym
            geometry=geom4.xyz
            {{multi;{cas}; wf,{wf};state,{state};start,2140.2;orbital,2243.2;}}
            {{mrci; {cas}; wf,{wf};state,{state};save,6003.2;noexc}}
            {{ci;trans,6000.2,6003.2;dm,8003.2}}


            !for -dp
            symmetry,nosym
            geometry=geom5.xyz
            {{multi;{cas}; wf,{wf};state,{state};start,2140.2;orbital,2244.2;}}
            {{mrci; {cas}; wf,{wf};state,{state};save,6004.2;noexc}}
            {{ci;trans,6000.2,6004.2;dm,8004.2}}


            '''.format( memory = self.memory,
                        basis = self.eInfo['basis'],
                        state=self.eInfo['state'],
                        wf =  self.eInfo['wf'],
                        cas = self.eInfo['cas']))


        nactTemp= ''

        for i,j in self.nactPairs:
            nactTemp+=textwrap.dedent(''' 
                !for taur     
                {{ddr,{dt},2140.2,2241.2,8001.2;state, {j}.1,{i}.1}}
                nacmepv=nacme

                {{ddr,-{dt},2140.2,2242.2,8002.2;state, {j}.1,{i}.1}}
                nacmemv=nacme

                nacmr = 0.5*(nacmepv+ nacmemv)

                !for taup
                {{ddr,{dp},2140.2,2243.2,8003.2;state, {j}.1,{i}.1}}
                nacmepv=nacme

                {{ddr,-{dp},2140.2,2244.2,8004.2;state, {j}.1,{i}.1}}
                nacmemv=nacme
                nacmp = 0.5*(nacmepv+ nacmemv)


                table, nacmr,nacmp
                save,ddrnact{i}{j}.res,new;

                
                '''.format(dt=self.dt,dp=self.dp,i=i,j=j))



        molproGridTemplate += nactTemp + '\n---\n'

        with open('grid.com','w') as f:
            f.write(molproGridTemplate)



    def makeGrid(self, ll):
        return np.arange(ll[0], ll[1]+ll[2], ll[2])



    def parseConfig(self, conFigFile):
        #parse configuration for running molpro
        scf = ConfigParser.SafeConfigParser()
        scf.read(conFigFile)

        molInfo =  dict(scf.items('molInfo'))
        self.scrdir = molInfo['scrdir']
        self.memory = molInfo['memory']

        self.eInfo = dict(scf.items('eInfo'))
        self.nInfo = dict(scf.items('nInfo'))
        gInfo = dict(scf.items('gInfo'))

        self.state = int(self.eInfo['state'])
        self.nactPairs = [[i,j] for i in range(1,self.state+1) for j in range(i+1,self.state+1)]

        #or should I just take as rho/theta/phi grid?
        self.rhoList = map(float, gInfo['rho'].split(','))
        self.thetaList = map(float, gInfo['theta'].split(','))
        self.phiList = map(float, gInfo['phi'].split(','))

        if self.nInfo['method']=='cpmcscf':
            if len(self.rhoList) == 1:
                raise Exception('Provide a rho grid for analytic calculation')
            self.rhoGrid = self.makeGrid(self.rhoList)
            self.phiGrid = self.makeGrid(self.phiList)
            mInfo = dict(scf.items('mInfo'))
            self.vModes = [int(i)-1 for i in mInfo['varying'].split(',')]
            self.createAnaTemplate()
            self.getTauThetaPhi = self.getTauThetaPhiAna


        elif self.nInfo['method']=='ddr':
            if len(self.rhoList) != 1:
                raise Exception('Give a fixed rho value for ddr calculation')
            self.rho = self.rhoList[0]
            self.thetaGrid = self.makeGrid(self.thetaList)
            self.phiGrid = self.makeGrid(self.phiList)
            self.dt = gInfo['dt']
            self.dp = gInfo['dp']
            self.createDdrTemplate()
            self.createGridGeom = self.createDdrGridGeom
            self.getTauThetaPhi = self.getTauThetaPhiDdr

        else:
            sys.exit('%s Not a proper option'%self.nInfo['method'])


    def parseResult(self, file):
        with open(file,"r") as f:
            dat = f.read().replace("D","E").strip().split("\n")[1:]
        dat = [map(float,i.strip().split()) for i in dat]
        return np.array(dat)



class Spectra(Base):

    def __init__(self, 
            conFigFile= 'molpro.config', 
            geomFile  = 'equiGeom.dat', 
            freqFile  = 'freq.dat', 
            wilsonFile= 'wilson.dat' ):



        self.parseData(geomFile, freqFile, wilsonFile)
        self.parseConfig(conFigFile)


    def parseData(self, geomFile, freqFile, wilsonFile):
        geomData = np.loadtxt(geomFile, 
            dtype={'names': ('names', 'mass', 'x','y','z'),'formats': ('S1', 'f8', 'f8','f8','f8')})
        atomNames= geomData['names']
        atomsMass= geomData['mass']
        equiGeom = np.array([geomData[i] for i in ['x','y','z']]).T

        freq = np.loadtxt(freqFile)
        wilson = np.loadtxt(wilsonFile)

        wilson = wilson.reshape(atomNames.shape[0], 3, freq.shape[0])
        freqInv = np.sqrt(self.hbar/(freq*self.cInvToTauInv))
        massInv = np.sqrt(1/atomsMass)
        wilFM = np.einsum('ijk,k,i->ijk',wilson,freqInv,massInv)
        self.atomNames= atomNames
        self.equiGeom = equiGeom
        self.wilFM    = wilFM



    def createGridGeom(self, rho, phi, outFile = 'geom.xyz'):
        
        nModes = self.wilFM.shape[2]
        qCord  = np.zeros(nModes)
        qCord[self.vModes[0]] = rho*self.cos(phi)
        qCord[self.vModes[1]] = rho*self.sin(phi)

        curGeom  = self.equiGeom+np.einsum('ijk,k->ij',self.wilFM, self.qCord)
        msg = 'for Rho = {}, Phi = {}'.format(rho, phi)
        nAtoms = len(self.atomNames)
        tmp = " {}\n".format(nAtoms)
        tmp+= "Geometry file created from ADT program. %s \n"%msg
        for i in range(nAtoms):
            tmp += "{},{},{},{}\n".format(self.atomNames[i], *curGeom[i])

        with open(outFile,"w") as f:
            f.write(tmp)


    def createDdrGridGeom(self, rho, phi):
        createGridGeom(atomNames, rho,  phi,  'geom1.xyz')
        createGridGeom(atomNames, rho+self.dr, phi, 'geom2.xyz')
        createGridGeom(atomNames, rho-self.dr, phi, 'geom3.xyz')
        createGridGeom(atomNames, rho, phi+self.dp,'geom4.xyz')
        createGridGeom(atomNames, rho, phi-self.dp,'geom5.xyz')
    

    def parseNact(self, i,j,rho,phi):
        file = 'ananac{}{}.res'.format(i,j)

        grads = self.angtobohr*parseResult(file)
        rotoMat = np.array([[self.cos(phi), self.sin(phi)], [-rho*self.sin(phi), rho*self.cos(phi)]])
        tau = np.einsum('ijk,ij,lk->l',self.wilFM[...,self.vModes], grads, rotoMat)
        return np.abs(tau)



    def getTauThetaPhiAna(self, rho, phi):    
        return np.stack((parseNact(i, j, rho, phi) 
                                    for i,j in self.nactPairs)).T


    def getTauThetaPhiDdr(self, *args):
        return np.stack((self.parseResult('ddrnact{}{}.res'.format(i,j)) 
                                    for i,j in self.nactPairs)).T



    def equiRun(self):
        # how to save this data??????
        nAtoms = len(self.atomNames)
        tmp = " {}\n".format(nAtoms)
        tmp+= "Geometry file created from ADT program. %s \n"
        for i in range(nAtoms):
            tmp += "{},{},{},{}\n".format(self.atomNames[i], *self.equiGeom[i])

        with open('geom.xyz',"w") as f:
            f.write(tmp)
        exitcode = subprocess.call(['molpro',"-d", self.scrdir ,'-W .','equi.com'])
        if exitcode==1: sys.exit('Molpro failed in equilibrium step')
        equiData = parseResult('equienr.res').flatten()


    def runMolpro(self):

        #initialise blank arrays to store result
        energyResult  = np.array([], dtype=np.float64).reshape(0,self.state+2) 
        nactRhoResult = np.array([], dtype=np.float64).reshape(0,nTau+2)
        nactPhiResult = np.array([], dtype=np.float64).reshape(0,nTau+2)


        #run molpro for the equilibrium step
        self.equiRun()

        for phi in self.phi_grid[:1]:
            shutil.copy('molpro_equi.wfu', 'molpro.wfu')   # copy the wave function from the equilibrium step
            for rho in self.rho_grid[:1]:

                self.createGridGeom(rho, phi) 
                shutil.copy('molpro.wfu',self.scrdir )

                exitcode = subprocess.call(['molpro',"-d", self.scrdir ,'-W .','grid.com'])
                if exitcode==0:
                    print 'Job successful  '+msg
                else : 
                    print 'Job unsuccessful'+msg 
                    continue


                enrData = parseResult('enr.res').flatten()
                tauRho, tauPhi = self.getTauThetaPhi(theta, phi)


                energyResult  = np.vstack((energyResult,  np.append([rho,phi],enrData)))
                nactRhoResult = np.vstack((nactRhoResult, np.append([rho,phi],tauRho)))
                nactPhiResult = np.vstack((nactPhiResult, np.append([rho,phi],tauPhi)))

        #or just save as npy format
        self.writeFile('energy.dat', energyResult)
        self.writeFile('tau_rho.dat', nactRhoResult)
        self.writeFile('tau_phi.dat', nactPhiResult)


    def writeFile(self, file, data):
        #along rho or along phi?
        for tp in self.rho_grid:
            np.savetxt( file, data[data[:,fc]==tp] ,delimiter="\t", fmt="%.8f")
            file.write("\n")



class ScatterFuncs(Base):

    #provide atom mass and names in atomFile.dat
    def __init__(self, 
            conFigFile= 'molpro.config',
            atomFile  = 'atomFile.dat'):

        self.parseData(atomFile)
        self.parseConfig(conFigFile)


    def AreaTriangle(self,a,b,c):
        """ area of a tringle with sides a,b,c """
        ps = (a+b+c)/2.0
        ar = ps*(ps-a)*(ps-b)*(ps-c)
        # negative area due to round off errors set to zero
        if ar < 0.0:
            ar = 0.0
        ar = math.sqrt(ar)
        return ar



    def to_jacobi(self,theta,phi):

        """ returns jacobi coordinates """
        m1, m2, m3 = self.atomMass

        M = m1 + m2 + m3
        mu = np.sqrt(m1*m2*m3/M)
        d1 = np.sqrt(m1*(m2+m3)/(mu*M))
        d2 = np.sqrt(m2*(m3+m1)/(mu*M))
        d3 = np.sqrt(m3*(m1+m2)/(mu*M))
        eps3 = 2 * np.atan(m2/mu)
        eps2 = 2 * np.atan(m3/mu)

        R1 = (1.0/np.sqrt(2.0))*self.rho*d3*np.sqrt(1.0+np.sin(theta)*np.cos(phi+eps3)) # F-H2 distance
        R2 = (1.0/np.sqrt(2.0))*self.rho*d1*np.sqrt(1.0+np.sin(theta)*np.cos(phi))      # H1-H2 distance
        R3 = (1.0/np.sqrt(2.0))*self.rho*d2*np.sqrt(1.0+np.sin(theta)*np.cos(phi-eps2)) # F-H1 distance

        if R1 < 1e-10:
            R1 = 0.0
        if R2 < 1e-10:
            R2 = 0.0
        if R3 < 1e-10:
            R3 = 0.0

        area = self.AreaTriangle(R1,R2,R3)
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
        rs = np.sqrt(np.dot(r,r))
        rc = np.sqrt(np.dot(R,R))
        if rc < 1e-10:
            rc = 0.0

        rtil = (R2*m2/(m2+m3),0.0)
        drtil = np.sqrt(np.dot(rtil,rtil))
        Areasmall = self.AreaTriangle(drtil,R1,rc)
        y = 4.0 * Areasmall
        x = drtil*drtil + rc*rc - R1*R1
        if np.fabs(x) < 1e-10:
            x = 0.0
        gamma = np.atan2(y,x)

        # we assert that gamma is always less than 90 degree always - our grid does not produce this for H2F
        # assert (gamma <= np.pi/2)

        return (rs, rc, gamma)



    def hyperToCart(self, theta, phi):
        rs, rc, gamma = to_jacobi(theta, phi)
        #! this gamma is in radian, so numpy sin cos can be used.
        p1 = [0, rc*np.cos(gamma),rc*np.sin(gamma)]
        p2 = [0, -rs/2.0 , 0.0 ]
        p3 = [0, rs/2.0  , 0.0 ]
        return np.array([p1,p2,p3])


    def parseData(self, atomFile):
        atomData = np.loadtxt(atomFile, 
            dtype={'names': ('names', 'mass'),'formats': ('S1', 'f8')})

        self.atomNames = atomData['names']
        self.atomMass  = atomData['mass']


    def createGridGeom(self, theta, phi, outFile = 'geom.xyz'):
        curGeom  = hyperToCart(theta, phi)
        msg = 'for Rho = {}, Theta={}, Phi = {}'.format(self.rho, theta, phi)
        nAtoms = len(self.atomNames)
        tmp = " {}\n".format(nAtoms)
        tmp+= "Geometry file created from ADT program. %s \n"%msg
        for i in range(nAtoms):
            tmp += "{},{},{},{}\n".format(self.atomNames[i], *curGeom[i])

        with open(outFile,"w") as f:
            f.write(tmp)


    def createDdrGridGeom(self, theta, phi):
        createGridGeom(atomNames, theta,    phi,    'geom1.xyz')
        createGridGeom(atomNames, theta+self.dt, phi,   'geom2.xyz')
        createGridGeom(atomNames, theta-self.dt, phi,   'geom3.xyz')
        createGridGeom(atomNames, theta,    phi+self.dp,'geom4.xyz')
        createGridGeom(atomNames, theta,    phi-self.dp,'geom5.xyz')
    

    def parseNact(self, i,j,rho,phi, gradTheta, gradPhi):
        file = 'ananac{}{}.res'.format(i,j)
        grads = self.angtobohr*parseResult(file)

        tauTheta = np.einsum('ij,ij', grads, gradTheta)
        tauPhi   = np.einsum('ij,ij', grads, gradPhi)
        return np.array([tauTheta, tauPhi])



    def getTauThetaPhiAna(self, theta, phi):
        # What will be this values
        dTheta = self.thetaList[2]/100
        dPhi = self.phiList[2]/100

        gradThetaPlus  = hyperToCart(theta+dTheta, phi)
        gradThetaMinus = hyperToCart(theta-dTheta, phi)
        gradTheta      = (gradThetaPlus - gradThetaMinus)/2*dTheta

        gradPhiPlus  = hyperToCart(theta+dTheta, phi)
        gradPhiMinus = hyperToCart(theta-dTheta, phi)
        gradPhi      = (gradPhiPlus - gradPhiMinus)/2*dPhi 
    
        return np.stack((parseNact(i,j,rho,phi, gradTheta, gradPhi) 
                                    for i,j in self.nactPairs)).T


    def getTauThetaPhiDdr(self, *args):
        return np.stack((self.parseResult('ddrnact{}{}.res'.format(i,j)) 
                                    for i,j in self.nactPairs)).T


    #include the theta =0 , phi = 0 step
    def equiRun(self):
        # run theta,phi 0,0 step as equilibrium step
        self.createGridGeom(self, 0.0, 0.0)

        exitcode = subprocess.call(['molpro',"-d", self.scrdir ,'-W .','equi.com'])
        if exitcode==1: sys.exit('Molpro failed in equilibrium step')
        equiData = parseResult('equienr.res').flatten()
        #what to do with this data?




    def runMolpro(self):

        #initialise blank arrays to store result
        energyResult  = np.array([], dtype=np.float64).reshape(0,self.state+2) 
        nactRhoResult = np.array([], dtype=np.float64).reshape(0,nTau+2)
        nactPhiResult = np.array([], dtype=np.float64).reshape(0,nTau+2)

        self.equiRun()

        for phi in self.phi_grid[:1]:
            shutil.copy('molpro_equi.wfu', 'molpro.wfu')
            for theta in self.theta_grid[:1]:

                self.createGridGeom(theta, phi) 
                shutil.copy('molpro.wfu',scrdir)

                exitcode = subprocess.call(['molpro',"-d", scrdir,'-W .','grid.com'])
                if exitcode==0:
                    print 'Job successful  '+msg
                else : 
                    print 'Job unsuccessful'+msg 
                    continue


                enrData = parseResult('enr.res').flatten()
                tauTheta, tauPhi = self.getTauThetaPhi(theta, phi)


                energyResult  = np.vstack((energyResult,  np.append([theta, phi],enrData)))
                nactThetaResult = np.vstack((nactThetaResult, np.append([theta, phi],tauTheta)))
                nactPhiResult = np.vstack((nactPhiResult, np.append([theta, phi],tauPhi)))

        self.writeFile('energy.dat', energyResult)
        self.writeFile('tau_theta.dat', nactThetaResult)
        self.writeFile('tau_phi.dat', nactPhiResult)


    def writeFile(self, file, data):
        #along rho or along phi?
        for tp in self.theta_grid:
            np.savetxt( file, data[data[:,fc]==tp] ,delimiter="\t", fmt="%.8f")
            file.write("\n")


if __name__ == "__main__":
    Spectra()

