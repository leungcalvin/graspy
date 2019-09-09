from grasp import *
# Todo: reproduce the following calculation:
#https://github.com/leungcalvin/vuletic/blob/3d0d463c11a5af4018d2a648f1db3602867eb892/172Yb/172Yb/script_1P1.out

testdir = '/home/calvin/grasp-wrapper/testdir2'
clistordering = ['1s','2s','2p','3s','3p','3d','4s','4p','4d','5s','5p','4f','5d','6s','6p','5f','6d','7p','8s']
initialize(workdir=testdir,clist=clistordering)
mr = Rcsfgenerate('Xe',
            ['4f(14,*)6s(1,*)6p(1,*)','4f(13,*)5d(1,*)6s(2,*)'],
            activeset=[6,6,5,4],
            jlower=0,jhigher=10,exc=0)

indices_6s6p = [[1],[1,2],[1],[],[],[]]
graspscript =[
        Rnucleus(Z=70,A=172,neutralMass=171.936378,I=0,NDM=0,NQM=0),
        Rangular(),
        Rwfnestimate(orbdict=
            {'*':'/home/calvin/Documents/vuletic/172Yb/172Yb/rwfn_6s2_6s6p.out'},
                     fallback='Thomas-Fermi'),
        Rmcdhf(indices_6s6p,orbs=['*'],specorbs=['*'],weightingmethod='Standard',runs=1000),
        Rcsfgenerate('Xe',
            ['4f(14,13)6s(1,*)6p(1,*)','4f(13,12)5d(1,*)6s(2,*)'],
            activeset=[7,6,5,4],
            jlower=0,jhigher=10,exc=1),
        Rcsfinteract('Dirac-Coulomb'),
        Rangular(),
        Rwfnestimate(orbdict={'*':os.path.join(testdir,'rwfn.out')},
            fallback='Thomas-Fermi'),
        Rmcdhf(indices_6s6p,orbs=['5d'],specorbs=['*'],weightingmethod='Standard',runs=1000),
        Rwfnestimate(orbdict={'*':os.path.join(testdir,'rwfn.out')},
            fallback='Thomas-Fermi'),
        Rmcdhf(indices_6s6p,orbs=['*'],specorbs=['*'],weightingmethod='Standard',runs=1000),
        Rwfnestimate(orbdict={'*':os.path.join(testdir,'rwfn.out')},
            fallback='Thomas-Fermi'),
        Rmcdhf(indices_6s6p,orbs=['7s'],specorbs=['*'],weightingmethod='Standard',runs=1000),
        Rsave('singlet-triplet'),
        Rci(calcname='singlet-triplet',
                includetransverse=True,
                modifyfreq=True,
                scalefactor='1.d-6',
                includevacpol = True,
                includenms= False,
                includesms = False,
                estselfenergy= True,
                largestn = 8,
                asfidx = indices_6s6p),
        ]

readoutscript = [
        Rmixextract('singlet-triplet',useCI=True,tolerance = 0.001,sort=True),
        JJtoLSJ('singlet-triplet',useCI=True,unique=True),
        Rlevels(calcname='singlet-triplet')
        ]
mr.execute(workdir = testdir, writeMR=True) # define the multireference
graspresult = [cmd.execute(workdir = testdir) for cmd in graspscript]
readoutresult = [cmd.execute(workdir = testdir) for cmd in readoutscript]
