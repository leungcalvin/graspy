from grasp import *
#This file expands the safronova calculation, starting with the results in '/home/calvin/graspy/safronova/ci'

clistordering = ['1s','2s','2p','3s','3p','3d','4s','4p','4d','5s','5p','4f','5d','6s','6p','5f','6d','7p','8s']

indices_mr_expanded = [[1],[1],
                       [1],[1,2,3],
                       [1,2],[1,2,3],
                       [1,2],[1,2],
                       [1],[1,2],
                       [ ],[1,2],
                           [1],
                           [ ],
                           [ ]]
dir_ci = '/home/calvin/graspy/safronova/ci'
dir_ci_expanded = '/home/calvin/graspy/safronova/ci_expanded'
initialize(workdir=dir_ci_expanded,clist=clistordering)
iccut_vals = Rmixaccumulate(calcname='safronova',useCI=True,truncate_eps=0.99).execute(workdir=dir_ci)
print(iccut_vals)
# generate 88765 expansion
singles_6s2_expanded = Rcsfgenerate('Xe',['4f(14,*)6s(2,*)'],activeset=[8,8,7,6,5],jlower=0,jhigher=2,exc=2)
singles_4f5d_expanded = Rcsfgenerate('Xe',['4f(14,*)5d(1,*)6s(1,*)'],activeset=[8,8,7,6,5],jlower=0,jhigher=10,exc=2)
singles_even_expanded = singles_6s2_expanded + singles_4f5d_expanded

singles_6s6p_expanded = Rcsfgenerate('Xe',['4f(14,*)6s(1,*)6p(1,*)',],activeset=[8,8,7,6,5],jlower=0,jhigher=6 ,exc=2)
singles_4f5d_expanded = Rcsfgenerate('Xe',['4f(13,*)5d(1,*)6s(2,*)',],activeset=[8,8,7,6,5],jlower=0,jhigher=16,exc=2)
singles_odd_expanded = singles_6s6p_expanded + singles_4f5d_expanded
singles_exp_expanded = singles_even_expanded + singles_odd_expanded
ci_run_expanded = [
         Rnucleus(Z=70,A=172,neutralMass=171.936378,I=0,NDM=0,NQM=0),
         singles_exp_expanded,
         Rcsfinteract('Dirac-Coulomb'),
         Rcsfzerofirst(small_exp = os.path.join(dir_ci,'rcsf.out'), big_exp = os.path.join(dir_ci_expanded,'rcsf.out')),
         Rangular(iccut=iccut_vals),

         Rwfnestimate(orbdict={'*' :os.path.join(dir_ci,'rwfn.out')},
                               fallback='Thomas-Fermi'),
         Rmcdhf(indices_mr_expanded,orbs=['5g*,6f*,7d*,8s,8p*'],specorbs=[''],weightingmethod='Standard',runs=1000),
         Rsave('safronova'),
         Rci(calcname='safronova',
                 includetransverse=True,
                 modifyfreq=True,
                 scalefactor='1.d-6',
                 includevacpol = True,
                 includenms= False,
                 includesms = False,
                 estselfenergy= True,
                 largestn = 8,
                 asfidx = indices_mr_expanded),

         ]

readout = [
        Rmixextract('safronova',useCI=True,tolerance = 0.001,sort=True),
        JJtoLSJ('safronova',useCI=True,unique=True),
        Rlevels(calcname='safronova')
        ]

#mr_6s2.execute(workdir = dir_6s2, writeMR=True) # define the multireference
#[cmd.execute(workdir = dir_6s2) for cmd in est6s2]
#mr_5d.execute(workdir = dir_5d, writeMR=True) # define the multireference
#[cmd.execute(workdir = dir_5d) for cmd in est5d]
#mr_6p.execute(workdir = dir_6p, writeMR=True) # define the multireference
#[cmd.execute(workdir = dir_6p) for cmd in est6p]

#multireference.execute(workdir=dir_ci,writeMR=True)
#[cmd.execute(workdir=dir_ci) for cmd in ci_run]

zeros_6s2 = Rcsfgenerate('Xe',['4f(14,*)6s(2,*)'],activeset=[6,6,6,5],jlower=0,jhigher=2,exc=0)
zeros_6s5d= Rcsfgenerate('Xe',['4f(14,*)5d(1,*)6s(1,*)'],activeset=[6,5,5,4],jlower=0,jhigher=8,exc=0)
zeros_even = zeros_6s2 + zeros_6s5d

zeros_6s6p= Rcsfgenerate('Xe',['4f(14,*)6s(1,*)6p(1,*)',],activeset=[6,6,4,4],jlower=0,jhigher=4,exc=0)
zeros_4f5d= Rcsfgenerate('Xe',['4f(13,*)5d(1,*)6s(2,*)',],activeset=[6,5,5,4],jlower=0,jhigher=14,exc=0)
zeros_odd = zeros_6s6p + zeros_4f5d

multireference = zeros_odd + zeros_even
multireference.execute(workdir=dir_ci_expanded,writeMR=True)

[cmd.execute(workdir=dir_ci_expanded) for cmd in ci_run_expanded]
[cmd.execute(workdir=dir_ci_expanded) for cmd in readout]
#        Rcsfinteract('Dirac-Coulomb'),
#        Rangular(),
#        Rwfnestimate(orbdict={'*':os.path.join(testdir,'rwfn.out')},
#            fallback='Thomas-Fermi'),
#        Rmcdhf(indices_6s6p,orbs=['5d'],specorbs=['*'],weightingmethod='Standard',runs=1000),
#        Rwfnestimate(orbdict={'*':os.path.join(testdir,'rwfn.out')},
#            fallback='Thomas-Fermi'),
#        Rmcdhf(indices_6s6p,orbs=['*'],specorbs=['*'],weightingmethod='Standard',runs=1000),
#        Rwfnestimate(orbdict={'*':os.path.join(testdir,'rwfn.out')},
#            fallback='Thomas-Fermi'),
#        Rmcdhf(indices_6s6p,orbs=['7s'],specorbs=['*'],weightingmethod='Standard',runs=1000),
#        ]

