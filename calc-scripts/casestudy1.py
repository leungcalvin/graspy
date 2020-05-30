from grasp import *
#
# This implements the GRASP 2018 script demonstration as found in the GRASP 2018 manual, using the GRASPy interface.
#
testdir = './casestudy1'
initialize(workdir=testdir)

# 1) Generate the odd multireference
odd_mr = Rcsfgenerate(core='None',ordering = 'Default',
            csflist=['1s(2,i)2s(2,i)2p(1,i)',
                     '1s(2,i)2p(3,i)'],
            activeset=[2,2],
            jlower=1,jhigher=3,exc=0)

odd_exp = Rcsfgenerate(core='None',ordering = 'Default',
            csflist=['1s(2,*)2s(2,*)2p(1,*)',
                     '1s(2,*)2p(3,*)'],
            activeset=[6,6,6,6,6,6],
            jlower=1,jhigher=3,exc=2)

odd_mr.execute(workdir = testdir, writeMR=True) #the execute() method performs the actual call to GRASP 2018.

odd_exp_split = Rcsfsplit(calcname = 'odd',nsplit=4,splitorbs = range(3,7))
odd_exp_split.execute(workdir = testdir)

# 2) Generate even CSF expansion split
even_2s2p2 = Rcsfgenerate(core='None',ordering = 'Default',
            csflist=['1s(2,i)2s(1,i)2p(2,i)',
                     ],
            activeset=[2,2],
            jlower=1,jhigher=5,exc=0)

even_2s2p2 = Rcsfgenerate(core='None',ordering = 'Default',
            csflist=['1s(2,*)2s(1,*)2p(2,*)',
                     ],
            activeset=[6,6,6,6,6,6],
            jlower=1,jhigher=5,exc=2)

even_2s2p2_split = Rcsfsplit(calcname = 'even',nsplit = 4, splitorbs = range(3,7))

assert(1 == 0)

# 2) Perform a SCF procedure to solve for the 1s,2s,2p orbitals.
MR_DHF =[
        Rnucleus(Z=42,A=96,neutralMass=96,I=1,NDM=1,NQM=1),
        Rangular(),
        Rwfnestimate(orbdict = None, fallback='Thomas-Fermi'),
        Rmcdhf([[1],[1],[1]],orbs = ['*'],specorbs = ['*'], runs = 100, weightingmethod = 'Standard'),
        Rsave('2s_2p_DF') ]
[cmd.execute(workdir = testdir) for cmd in MR_DHF]

# 3) Generate a CAS expansion from the 2S configuration.
CAS_2S_exp = Rcsfgenerate(core='None',ordering = 'Default',
            csflist=['1s(2,*)2s(1,*)'],
            activeset=[3,3,3],
            jlower=1,jhigher=1,exc=3)
CAS_2S_exp.execute(workdir = testdir)

# 4) Solve for the n=3 correlation orbitals, using orbitals generated from 2s_2p_DF.w above.
indices_2S = [[1]]
CAS_2S = [
        Rangular(),
        Rwfnestimate(orbdict = {'*':'2s_2p_DF.w'},fallback = 'Thomas-Fermi'),
        Rmcdhf(asfidx = indices_2S,
            orbs = ['3*'],
            specorbs = [' '],
            runs = 100, weightingmethod = 'Standard'),
        Rsave('2s_3')
        ]
[cmd.execute(workdir = testdir) for cmd in CAS_2S]

# 5) Perform CI on the 2S expansion.
CI_2S = [
        Rci(calcname='2s_3',
                includetransverse=True,
                modifyfreq=True,
                scalefactor='1.d-6',
                includevacpol = True,
                includenms= False,
                includesms = False,
                estselfenergy= True,
                largestn = 3,
                asfidx = indices_2S),
        JJtoLSJ(calcname= '2s_3',useCI = True, unique = True),
        ]
[cmd.execute(workdir = testdir) for cmd in CI_2S]

# 6) Generate a CAS expansion from the 2P configuration.
CAS_2P_exp = Rcsfgenerate(core='None',ordering = 'Default',
            csflist=['1s(2,*)2p(1,*)'],
            activeset=[3,3,3],
            jlower=1,jhigher=3,exc=3)

CAS_2P_exp.execute(workdir = testdir, writeMR = False)

# 7) Solve for the n=3 correlation orbitals, using orbitals generated from 2s_2p_DF.w above.
indices_2P = [[1],[1]]
CAS_2P = [
        Rangular(),
        Rwfnestimate(orbdict = {'*':'2s_2p_DF.w'},fallback = 'Thomas-Fermi'),
        Rmcdhf(asfidx = indices_2P,orbs = ['3*'],specorbs = [' '],runs = 100, weightingmethod = 'Standard'),
        Rsave('2p_3')
        ]
[cmd.execute(workdir = testdir) for cmd in CAS_2P]

# 8) Perform CI on the 2P expansion.
CI_2P = [
        Rci(calcname='2p_3',
                includetransverse=True,
                modifyfreq=True,
                scalefactor='1.d-6',
                includevacpol = True,
                includenms= False,
                includesms = False,
                estselfenergy= True,
                largestn = 3,
                asfidx = indices_2P),
        JJtoLSJ(calcname= '2p_3',useCI = True, unique = True),
        ]
[cmd.execute(workdir = testdir) for cmd in CI_2P]

# 9) Calculate transitions.
transitions_2P = [
        Rhfs(calcname = '2p_3',useCI=True),
        Rbiotransform(useCI=True,calcname_initial = '2s_3',calcname_final = '2p_3', transform_all = True),
        Rtransition(useCI=True,calcname_initial = '2s_3',calcname_final = '2p_3',transition_spec = ['E1'])]

[cmd.execute(workdir = testdir) for cmd in transitions_2P]
