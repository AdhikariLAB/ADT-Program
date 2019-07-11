
# wf now splitted into electron, spin and charge numbers and the wf card will be built inside the scripts
#? the symmetry information is inside the sysInfo section ????
# so eInfo now has this mandatory keyworsd

# 1. method    ----> same as before
# 2. basis     ----> same as before
# 3. cas       ----> occ,2;closed,2         no symmetry
#                    occ,2,3;closed,2,3     Cs
#                    occ,2,3,4,5;...........C2v so on...
# 4. electron  ----> number of electron, one number
# 5. spin      ----> one number
# 6. charge    ----> one number
# 7. state     ----> n numbers seperated by comma, n= number of IREPS 
#                    0 is mandatory for no state calculation for a particular state




###### An these optional keywords
# 1. scale 
# 2. uhf extra     ---> what will be done with this when wf is provided through this ?????
# 3. multi extra 
# 4. mrci extra




# def validity(self, cas, subs, nIreps):
#     ss = re.search('[0-9a-zA-Z,;](%s[0-9,]+;?)'%subs, cas).group(1)
#     assert len(re.findall('(\d+)', ss))==nIreps, "{} number of {} required for {} symmetry".format(nIreps, subs, symmetry)




def createSymAnaTempl(self):


    # number of IREPS    <<<=====
    # nIREPs = 
    # but occ closed and core is passed through the cas option so no manual prep

    # core can have one or more number due to number of ireps
    

    #* sanity checks for cas, check number of things like occ,closed, core is consistent with IREPs
    for substring in self.eInfo['cas'].split(';'):
        for item in ['occ', 'closed', 'core']:
            if item in substring:
                assert len(re.findall('(\d+)', substring)) == nIREPs, (
                    "{a} number of {b} required for {a} IREPs of {c} symmetry".format(
                                                a = nIREPs,
                                                b = item,
                                                c = self.symmetry
                                                )
                    )


    #* sanity checks for state and split state into list of integers
    nStateList = [int(i) for i in re.findall('(\d+)', self.eInfo['state'])]
    assert nStateList==nIREPs, (
            "{a} number of states required for {a} IREPs of {c} symmetry".format(
                                        a = nIREPs,
                                        c = self.symmetry
                                        )
            )


    cas = re.sub('[0-9a-zA-Z,;](core[0-9,]+;?)','', self.eInfo['cas'])


    irepCard = ''
    #looping through the statelist as its already checked through in above
    for nIrep, state in enumerate(nStateList, start=1):
        if(state):
            irepCard += 'wf,{elec},{irep},{spin},{chrg};state,{stat};'.format(
                            elec = self.eInfo['nElectron'],
                            irep = nIrep,
                            spin = self.eInfo['nSpin'],
                            chrg = self.eInfo['nCharge'],
                            stat = state
            )



    energyLine ='''{{mcscf;{cas};{irepCard};start,2140.2; orbital,2140.2;{extra}}}
                '''.format(
                        cas = cas   # <<<----- the cas after stripping the core 
                        irepCard =  irepCard,
                        extra= self.eInfo['multi_extra']
                        )



    if self.eInfo['method'] == 'mrci':
        energyLine+='''{{mrci;{cas}; {irepCard};{extra}}}
                    '''.format(
                            cas=self.eInfo['cas'],
                            irepCard = irepCard,
                            extra = self.eInfo['mrci_extra']
                            )




    molproTemplate=textwrap.dedent('''
        ***, Molpro template created from ADT program for analytical job.
        memory,{memory}
        file,2,molpro.wfu;

        basis={basis};

        symmetry,{sym}

        geometry=geom.xyz
        {enrl}

        show, energy
        table, energy
        save,enr.res,new
        {{table,____; noprint,heading}}
        '''.format(memory = self.memory,
                    basis = self.eInfo['basis'],
                    sym = self.sysInfo['symmetry']
                    enrl = energyLine))
        

    ## NACTs
    # nact pairs are now 3 dimensional list. i.e list containing the list of pairs
    # take it as states variable
    self.nactPairsList =[ [[i, j] for j in range(2, state+1) for i in range(1, j)] for state in nStateList]

    self.nTau = [len(item) for item in self.nactPairsList]

    nactText= "\n\nbasis={}\n".format(self.nInfo['basis'])


    getChunk = lambda seq: [seq[i:i+5] for i in range(0, len(seq), 5)]

    for nIrep, nactPairs in enumerate(self.nactPairsList, start=1):
        # check for 0 state
        # nactPairs contain the list of nact pairs for the irep --> nIrep
        # now from those pairs, have to calculate at max 5 pairs at a time
        for pairs in getChunk(nactPairs):
            nactText += "\n{{mcscf;{cas}; {irepCard}; start,2140.2;{extra};\n".format(
                    cas=cas,
                    irepCard = irepCard
                    extra= self.eInfo['multi_extra'])
            
            forceText = ''
            for count,pair in enumerate(pairs start=1):
                nactText += "cpmcscf,nacm,{f}.{irep},{s}.{irep},record=51{n:02}.1,{extra};\n".format(
                            f=f, s=s, 
                            n=count, 
                            irep = nIrep,
                            extra=self.nInfo['nact_extra'])

                #! CAUTION  Name of the res files         <<<< =============
                forceText +=textwrap.dedent("""
                            force;nacm,51{n:02}.1;varsav
                            table,gradx,grady,gradz;
                            save,ananac{f}{s}.res,new;
                            {{table,____;  noprint,heading}}

                            """.format(f=f, s=s, n=count))
            nactText += "}\n\n" + forceText




    molproTemplate += nactText + '\n---\n'

    with open('grid.com','w') as f:
        f.write(molproTemplate)






def createDdrTemplate(self):
    ''''Creates the molpro template files ddr job'''


    cas = re.sub('[0-9a-zA-Z,;](core,\d+;*)', '', self.eInfo['cas'])

    #mrci has to be done for ddr nact calculation
    energyLine = """
        {{mcscf;{cas1}; wf,{wf};state,{state};start,2140.2; orbital,2140.2;{extra1}}}
        {{mrci; {cas}; wf,{wf};state,{state};save,6000.2;dm,8000.2;{extra2}}}
        """.format(state   =self.eInfo['state'],
                    wf     = self.eInfo['wf'],
                    cas    = self.eInfo['cas'],
                    cas1   = cas,
                    extra1 = self.eInfo['multi_extra'],
                    extra2 = self.eInfo['mrci_extra'])

    # d1 is increment in rho, theta or q depending on calculation type
    # and d2 is always phi
    molproTemplate=textwrap.dedent('''
        ***, Molpro template created from ADT program for analytical job for analytical job.
        memory,{memory}
        file,2,molpro.wfu;

        basis={basis};

        symmetry,nosym

        geomtyp=xyz
        geometry=geom1.xyz

        {enrl}

        show, energy
        table, energy
        save,enr.res,new
        {{table,____; noprint,heading}}

        !for +d1
        symmetry,nosym
        geometry=geom2.xyz
        {{multi;{cas1} ;wf,{wf};state,{state};start,2140.2;orbital,2241.2;{extra1}}}
        {{mrci; {cas}; wf,{wf};state,{state};save,6001.2;{extra2}}}
        {{ci;trans,6000.2,6001.2;dm,8001.2;{extra3}}}


        !for -d1
        symmetry,nosym
        geometry=geom3.xyz
        {{multi;{cas1}; wf,{wf};state,{state};start,2140.2;orbital,2242.2;{extra1}}}
        {{mrci; {cas}; wf,{wf};state,{state};save,6002.2;{extra2}}}
        {{ci;trans,6000.2,6002.2;dm,8002.2;{extra3}}}



        !for +d2
        symmetry,nosym
        geometry=geom4.xyz
        {{multi;{cas1}; wf,{wf};state,{state};start,2140.2;orbital,2243.2;{extra1}}}
        {{mrci; {cas}; wf,{wf};state,{state};save,6003.2;{extra2}}}
        {{ci;trans,6000.2,6003.2;dm,8003.2;{extra3}}}


        !for -d2
        symmetry,nosym
        geometry=geom5.xyz
        {{multi;{cas1}; wf,{wf};state,{state};start,2140.2;orbital,2244.2;{extra1}}}
        {{mrci; {cas}; wf,{wf};state,{state};save,6004.2;{extra2}}}
        {{ci;trans,6000.2,6004.2;dm,8004.2;{extra3}}}


        '''.format( memory = self.memory,
                    basis  = self.eInfo['basis'],
                    enrl   = energyLine,
                    state  = self.eInfo['state'],
                    wf     = self.eInfo['wf'],
                    cas    = self.eInfo['cas'],
                    cas1   = cas,
                    extra1 = self.eInfo['multi_extra'],
                    extra2 = self.eInfo['mrci_extra'],
                    extra3 = self.nInfo['nact_extra']))


    nactTemp= ''

    for i,j in self.nactPairs:

        #implementing three point cebtral difference
        nactTemp+=textwrap.dedent(''' 
            !for taur     
            {{ddr, 2*{d1}
            orbital,2140.2,2141.2,2142.2;
            density,8000.2,8001.2,8002.2;
            state, {j}.1,{i}.1
            }}
            nacmr = nacme

            !for taup
            {{ddr, 2*{d2}
            orbital,2140.2,2143.2,2144.2;
            density,8000.2,8003.2,8004.2;
            state, {j}.1,{i}.1
            }}
            nacmp = nacme

            table, nacmr,nacmp
            save,ddrnact{i}{j}.res,new;

            '''.format(d1=self.d1,d2=self.d2,i=i,j=j))


    molproTemplate += nactTemp + '\n---\n'

    with open('grid.com','w') as f:
        f.write(molproTemplate)

    molproInitTemplate = textwrap.dedent('''
        ***, Molpro template created from ADT program for ddr job.
        memory,{memory}
        file,2,molpro_init.wfu,new;

        basis={basis};

        symmetry,nosym

        geomtyp=xyz
        geometry=geom.xyz

        {{hf}}
        {enrl}

        show, energy
        table, energy
        save,equienr.res,new
        {{table,____; noprint,heading}}

        '''.format( memory = self.memory,
                    basis = self.eInfo['basis'],
                    enrl=energyLine,
                    hf = self.eInfo['hf']))



    with open('init.com', 'w') as f:
        f.write(molproInitTemplate)




    molproScaleTemplate = textwrap.dedent('''
        ***, Molpro template created from ADT program for ddr job.
        memory,{memory}
        file,2,molpro_init.wfu,new;

        basis={basis};

        symmetry,nosym

        geomtyp=xyz
        geometry=geom.xyz

        {enrl}

        show, energy
        table, energy
        save,scale.res,new
        {{table,____; noprint,heading}}

        '''.format( memory = self.memory,
                    basis = self.eInfo['basis'],
                    enrl=energyLine))

    with open('scale.com', 'w') as f:
        f.write(molproScaleTemplate)

