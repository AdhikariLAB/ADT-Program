
def createTemplate(eInfo, nInfo):
    molproEquiTemplate=textwrap.dedent('''
        ***, Molpro template created from ADT program for "WhatEver"
        file,2,molpro.wfu,new;

        basis=6-31G**;

        symmetry,nosym


        geometry=geom.xyz

        {{uhf;orbital,2200.2}}
        {{{method};{cas}; wf,{wf};state,{state};orbital,2140.2;}}

        show, energy
        table, energy
        save,equienr.res,new
        {{table,____; noprint,heading}}

        ---

        '''.format(method=eInfo['method'],
                    state=eInfo['state'],
                    wf =  eInfo['wf'],
                    cas = eInfo['cas']))

    with open('equi.com' ,'w') as f:
        f.write(molproEquiTemplate)



    molproGridTemplate=textwrap.dedent('''
        ***, Molpro template created from ADT program for "WhatEver"
        memory,100,m
        file,2,molpro.wfu;

        basis=6-31G**;

        symmetry,nosym


        geometry=geom.xyz

        {{{method};{cas}; wf,{wf};state,{state};start,2140.2; orbital,2140.2;}}

        show, energy
        table, energy
        save,enr.res,new
        {{table,____; noprint,heading}}
        '''.format(method=eInfo['method'],
                    state=eInfo['state'],
                    wf =  eInfo['wf'],
                    cas = eInfo['cas']))


    state = int(eInfo['state'])
    ########################
    #   Following is only for analytical NACT, 
    if nInfo['method'] == 'analytical': 
        nInfo['nmethod'] = 'cpmcscf'

    nactPairs = [[i,j] for i in range(1,state+1) for j in range(i+1,state+1)]


    nactTemp= "\n\nbasis={}\n".format(eInfo['basis'])

    for ind in range(0,len(nactPairs),5):
        nactTemp += "\n{{{method};{cas}; wf,{wf};state,{state}; start,2140.2;\n".format(
                            method=eInfo['method'],
                            state=eInfo['state'],
                            wf =  eInfo['wf'],
                            cas = eInfo['cas'])

        pairChunk = nactPairs[ind:ind+5]
        forceTemp = ''
        for i,j in pairChunk:
            nactTemp += "{nmethod},nacm,{f}.1,{s}.1,record=51{f}{s}.1;\n".format(
                        nmethod = nInfo['nmethod'],
                        f=i,s=j)
            forceTemp +=textwrap.dedent("""
            force;nacm,51{f}{s}.1;varsav
            table,gradx,grady,gradz;
            save,ananac{f}{s}.res,new;
            {{table,____;  noprint,heading}}

            """.format(f=i,s=j))
        nactTemp += "}\n" + forceTemp




    molproGridTemplate += nactTemp + '\n---\n'

    with open('grid.com','w') as f:
        f.write(molproGridTemplate)

    return state, nactPairs






def createGeometry(atomNames, equiGeom,wilFM,vModes, rho, phi):
    nModes = wilFM.shape[2]
    qCord  = np.zeros(nModes)
    qCord[vModes[0]] = rho*np.cos(phi)
    qCord[vModes[1]] = rho*np.sin(phi)

    curGeom  = equiGeom+np.einsum('ijk,k->ij',wilFM, qCord)
    
    nAtoms = len(atomNames)
    tmp = " {}\n".format(nAtoms)
    tmp+= "Geometry file created from ADT program. %s \n"%msg
    for i in range(nAtoms):
        tmp += "{},{},{},{}\n".format(atomNames[i], *geomData[i])

    with open('geom.xyz',"w") as f:
        f.write(tmp)