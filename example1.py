from grasp import *

# This implements the GRASP 2018 calculation for 1s2 2s 2S and 1s2 2p 2P in Li I, as found in the GRASP 2018 manual.

testdir = '/home/calvin/graspy/example1'
initialize(workdir=testdir)

# rcsfgenerate generates list of CSFs for 2S and 2P with 3 CSFs.

ref1 = Rcsfgenerate(core='None',ordering = 'Default',
            csflist=['1s(2,i)2s(1,i)'],
            activeset=[2],
            jlower=1,jhigher=1,exc=0)
ref2 = Rcsfgenerate(core='None',ordering = 'Default',
            csflist=['1s(2,i)2p(1,i)'],
            activeset = [1,2],
            jlower=1,jhigher=3,exc=0)
mr = ref1 + ref2
mr.execute(workdir = testdir, writeMR=True) # define the multireference
MR_DHF =[
        Rnucleus(Z=3,A=7,neutralMass=6.941,I=1.5,NDM=3.2564268,NQM=-0.04),
        Rangular(),
        Rwfnestimate(orbdict = None, fallback='Thomas-Fermi'),
        Rmcdhf([[1],[1],[1]],orbs = ['*'],specorbs = ['*'], runs = 100, weightingmethod = 'Standard'),
        Rsave('2s_2p_DF')
        ]
[cmd.execute(workdir = testdir) for cmd in MR_DHF]

CAS_2S_exp = Rcsfgenerate(core='None',ordering = 'Default',
            csflist=['1s(2,*)2s(1,*)'],
            activeset=[3,3,3],
            jlower=1,jhigher=1,exc=3)
CAS_2S_exp.execute(workdir = testdir)

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
CAS_2P_exp = Rcsfgenerate(core='None',ordering = 'Default',
            csflist=['1s(2,*)2p(1,*)'],
            activeset=[3,3,3],
            jlower=1,jhigher=3,exc=3)

CAS_2P_exp.execute(workdir = testdir, writeMR = False)

indices_2P = [[1],[1]]
CAS_2P = [
        Rangular(),
        Rwfnestimate(orbdict = {'*':'2s_2p_DF.w'},fallback = 'Thomas-Fermi'),
        Rmcdhf(asfidx = indices_2P,orbs = ['3*'],specorbs = [' '],runs = 100, weightingmethod = 'Standard'),
        Rsave('2p_3')
        ]
[cmd.execute(workdir = testdir) for cmd in CAS_2P]

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
