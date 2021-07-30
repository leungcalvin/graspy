from graspy.grasp import *
#ex4_workdir = '/home/cleu/example4_mpi'
ex4_workdir = 'home/arjunb/calc-outputs/example4_mpi'
initialize(workdir=ex4_workdir)

Rnucleus(Z=26,A=56,neutralMass=55.845,I=1,NDM=1,NQM=1).execute(workdir=ex4_workdir)
Rcsfgenerate(core='He',ordering='Default',
        csflist=['2s(2,i)2p(6,i)3s(2,i)',
                 '2s(2,i)2p(6,i)3p(2,i)',
                 '2s(2,i)2p(6,i)3s(1,i)3d(1,i)'],
        active_set=[3,3,3],
        jlower=0,jhigher=6,exc=0,write_csf='evenmr.inp').execute(workdir=ex4_workdir)
Rangular().execute(workdir=ex4_workdir)
Rwfnestimate(orbdict = None, fallback='Thomas-Fermi').execute(workdir=ex4_workdir)
Rmcdhf([[1,2,3],[1,2],[1,2,3,4],[1]],orbs = ['*'],specorbs = ['*'], weighting_method = 'Standard', runs = 100).execute(workdir = ex4_workdir)
Rsave('evenmr').execute(workdir=ex4_workdir)


Rcsfgenerate(core='He',ordering='Default',
        csflist=['2s(2,1)2p(6,i)3s(2,*)',
        '2s(2,i)2p(6,5)3p(2,*)',
        '2s(2,1)2p(6,i)3p(2,*)',
        '2s(2,i)2p(6,5)3s(1,*)3d(1,*)',
        '2s(2,1)2p(6,i)3s(1,*)3d(1,*)'],
        active_set=[4,4,4,4],
        jlower=0,jhigher=6,exc=2,write_csf = 'even4.inp').execute(workdir=ex4_workdir)
Rangular().execute_mpi(workdir=ex4_workdir,nproc=4)
Rwfnestimate(orbdict = {'*':'rwfn.out'},fallback = 'Thomas-Fermi').execute(workdir=ex4_workdir)

indices_even = [[1,2,3],[1,2],[1,2,3,4],[1]] # states we care about
Rmcdhf(indices_even,orbs=['4*'],specorbs=[''],weighting_method='Standard',runs=100).execute_mpi(workdir=ex4_workdir,nproc=4)
#uncommented rmcdhf, might need to change back to commented
Rsave('even4').execute(workdir=ex4_workdir)

Rci(calc_name='even4',
        include_transverse = True,
        modify_freq = True,
        scale_factor = '1.d-6',
        include_vacpol = True,
        include_nms = False,
        include_sms = False,
        est_self_energy = True,
        largest_n = 4,
        asfidx = indices_even).execute_mpi(workdir = ex4_workdir)
JJtoLSJ(calc_name='even4',use_ci = True, unique = True).execute(workdir=ex4_workdir)

Rcsfgenerate(core='He',ordering='Default',
        csflist=['2s(2,1)2p(6,i)3s(1,i)3p(1,i)',
                 '2s(2,i)2p(6,i)3p(1,i)3d(1,i)'],
        active_set=[3,3,3],
        jlower=0,jhigher=8,exc=0,write_csf = 'oddmr.inp').execute(workdir=ex4_workdir)
Rangular().execute(workdir=ex4_workdir)
Rwfnestimate(orbdict = None, fallback='Thomas-Fermi').execute(workdir=ex4_workdir)
indices_odd = [[1,2],[1,2,3,4,5],[1,2,3,4,5],[1,2,3],[1]]
Rmcdhf(indices_odd,orbs = ['*'],specorbs = ['*'], weighting_method = 'Standard', runs = 100).execute(workdir = ex4_workdir)
Rsave('oddmr').execute(workdir=ex4_workdir)
Rcsfgenerate(core='He',ordering='Default',
        csflist=['2s(2,i)2p(6,5)3s(1,*)3p(1,*)',
        '2s(2,1)2p(6,i)3s(1,*)3p(1,*)',
        '2s(2,i)2p(6,5)3p(1,*)3d(1,*)',
        '2s(2,1)2p(6,i)3p(1,*)3d(1,*)'],
        active_set=[4,4,4,4],
        jlower=0,jhigher=8,exc=2,write_csf = 'odd4.inp').execute(workdir=ex4_workdir)
Rangular().execute_mpi(workdir=ex4_workdir,nproc=4)
Rwfnestimate(orbdict = {'*':'rwfn.out'},fallback = 'Thomas-Fermi').execute(workdir=ex4_workdir)

Rmcdhf(indices_odd,orbs=['4*'],specorbs=[''],weighting_method='Standard',runs=100).execute_mpi(workdir=ex4_workdir,nproc=4)
Rsave('odd4').execute(workdir=ex4_workdir)
Rci(calc_name='odd4',
        include_transverse = True,
        modify_freq = True,
        scale_factor = '1.d-6',
        include_vacpol = True,
        include_nms = False,
        include_sms = False,
        est_self_energy = True,
        largest_n = 4,
        asfidx = indices_odd).execute(workdir = ex4_workdir)
JJtoLSJ(calc_name='odd4',use_ci = True, unique = True).execute(workdir=ex4_workdir)
