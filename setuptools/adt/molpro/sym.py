
# wf now splitted into electron, spin and charge numbers and the wf card will be built inside the scripts
# the symmetry information is inside the sysInfo section ????
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




def validity(self, cas, subs, nIreps):
    ss = re.search('[0-9a-zA-Z,;](%s[0-9,]+;?)'%subs, cas).group(1)
    assert len(re.findall('(\d+)', ss))==nIreps, "{} number of {} required for {} symmetry".format(nIreps, subs, symmetry)




def createSymAnaTempl(self):


    


    # number of IREPS    <<<=====
    # nIREPs = 
    # but occ closed and core is passed through the cas option so no manual prep

    # core can have one or more number due to number of ireps
    cas = re.sub('[0-9a-zA-Z,;](core[0-9,]+;?)','', self.eInfo['cas'])

    #* sanity checks for cas, check number of things like occ,closed, core is consistent with IREPs




    #* sanity checks for state and split state into list of integers



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