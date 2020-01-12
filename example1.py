from grasp import *
#
# This implements the GRASP 2018 calculation for 1s2 2s 2S and 1s2 2p 2P in Li I, as found in the GRASP 2018 manual, using the GRASPy interface.
#
testdir = './example1-test'
initialize(workdir=testdir)

# 1) Generate a multireference consisting of the 1s2 2s and 1s2 2p configurations
ref_2s = Rcsfgenerate(core='None',ordering = 'Default',
            csflist=['1s(2,i)2s(1,i)'],
            activeset=[2],
            jlower=1,jhigher=1,exc=0)
ref_2p = Rcsfgenerate(core='None',ordering = 'Default',
            csflist=['1s(2,i)2p(1,i)'],
            activeset = [1,2],
            jlower=1,jhigher=3,exc=0)

# The three states which make up the multireference can be added together with arithmetic.
mr = ref_2s + ref_2p
mr.execute(workdir = testdir, writeMR=True) #the execute() method performs the actual call to GRASP 2018.

# 2) Perform a SCF procedure to solve for the 1s,2s,2p orbitals.
MR_DHF =[
        Rnucleus(Z=3,A=7,neutralMass=6.941,I=1.5,NDM=3.2564268,NQM=-0.04),
        Rangular(),
        Rwfnestimate(orbdict = None, fallback='Thomas-Fermi'),
        Rmcdhf([[1],[1],[1]],orbs = ['*'],specorbs = ['*'], runs = 100, weightingmethod = 'Standard'),
        Rsave('2s_2p_DF')
        ]
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
