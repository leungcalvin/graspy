from graspy.grasp import *
import os
#This file implements the GRASP calculation documented here:

def full_calculation(calc_dir,active_set,exc,n_open):
    dir_6s = os.path.join(calc_dir,'est_6s2')
    dir_6p = os.path.join(calc_dir,'est_6p')
    dir_5d = os.path.join(calc_dir,'est_5d')
    dir_ci = os.path.join(calc_dir,f'ci_{n_open}_{"".join(map(str,active_set))}_x{exc}')

    clistordering = ['1s','2s','2p','3s','3p','3d','4s','4p','4d','5s','5p','4f','5d','6s','6p','5f','6d','7p','8s']
    initialize(workdir=dir_6s,clist=clistordering)
    mr_6s2 = Rcsfgenerate('Xe',
                ['4f(14,*)6s(2,*)'],
                active_set=[6,5,4,4],
                jlower=0,jhigher=0,exc=0,write_csf='rcsfmr.inp')
    indices_6s2 = [[1]]
    est6s2 =[
            Rnucleus(Z=70,A=172,neutralMass=171.936378,I=0,NDM=0,NQM=0),
            Rangular(),
            Rwfnestimate(orbdict=
                {'*':'../../../calc-scripts/cores/yb_6s2master.w'},
                         fallback='Thomas-Fermi'),
            Rmcdhf(indices_6s2,orbs=['6s'],specorbs=['*'],weighting_method='Standard',runs=1000),
            Rwfnestimate(orbdict=
                {'*':os.path.join(dir_6s,'rwfn.out')},
                         fallback='Thomas-Fermi'),
            Rmcdhf(indices_6s2,orbs=['4f*'],specorbs=['*'],weighting_method='Standard',runs=1000),
            Rwfnestimate(orbdict=
                {'*':os.path.join(dir_6s,'rwfn.out')},
                         fallback='Thomas-Fermi'),
            Rmcdhf(indices_6s2,orbs=['*'],specorbs=['*'],weighting_method='Standard',runs=1000),
            ]

#
# Estimate 6p orbitals using 4f13 6s2 6p configurations
#
    initialize(workdir=dir_6p,clist=clistordering)
    mr_6p = Rcsfgenerate('Xe',
                ['4f(13,*)6s(2,*)6p(1,*)'],
                active_set=[6,6,4,4],
                jlower=0,jhigher=10,exc=0,write_csf = 'rcsfmr.inp')
    indices_6p = [[1],[1,2,3],[1,2,3,4],[1,2,3],[1]]
    est6p =[
            Rnucleus(Z=70,A=172,neutralMass=171.936378,I=0,NDM=0,NQM=0),
            Rangular(),
            Rwfnestimate(orbdict=
                {'*':os.path.join(dir_6s,'rwfn.out')},
                         fallback='Thomas-Fermi'),
            Rmcdhf(indices_6p,orbs=['6p*'],specorbs=['*'],weighting_method='Standard',runs=1000)
            ]


#
# Estimate 5d orbitals using 4f13 6s2 5d configurations
#
    def estimate_5d(which_csf):
        if which_csf == 'odd':
            mr_5d = Rcsfgenerate('Xe',
                        ['4f(13,*)5d(1,*)6s(2,*)'],
                        active_set=[6,5,5,4],
                        jlower=0,jhigher=14,exc=0,write_csf = 'rcsfmr.inp')
            indices_5d = [[1],[1,2,3],[1,2,3,4],[1,2,3,4],[1,2,3,4],[1,2,3],[1]]
        if which_csf == 'even':
            mr_5d = Rcsfgenerate('Xe',
                        ['4f(14,*)5d(1,*)6s(1,*)'],
                        active_set=[6,5,5,4],
                        jlower=2,jhigher=6,exc=0,write_csf = 'rcsfmr.inp')
            indices_5d = [[1],[1,2],[1]] # just the singlet and triplet
        est5d =[
                mr_5d,
                Rnucleus(Z=70,A=172,neutralMass=171.936378,I=0,NDM=0,NQM=0),
                Rangular(),
                Rwfnestimate(orbdict=
                    {'*':os.path.join(dir_6s,'rwfn.out')},
                             fallback='Thomas-Fermi'),
                Rmcdhf(indices_5d,orbs=['5d*'],specorbs=['*'],weighting_method='Standard',runs=1000,integration_method = 3)
                ]
        return est5d

   ##################
   #
   #
   # DHF Calculations
   #
   ##################
    mr_6s2.execute(workdir = dir_6s) # define the multireference
    [cmd.execute(workdir = dir_6s) for cmd in est6s2]
    # input('6s OK?')
    mr_6p.execute(workdir = dir_6p) # define the multireference
    [cmd.execute(workdir = dir_6p) for cmd in est6p]
    # input('6p OK?')
    initialize(workdir=dir_5d,clist=clistordering)
    [cmd.execute(workdir = dir_5d) for cmd in estimate_5d('odd')]
    # input('5d OK?')

   ##################
   #
   #
   # CI Calculations
   #
   ##################

    AR_PREFIX = '3d(10,c)4s(2,*)4p(6,*)4d(10,*)5s(2,*)5p(6,*)'
    def gen_multireference():
        zeros_6s2 = Rcsfgenerate('Ar',[AR_PREFIX + '4f(14,*)6s(2,*)'],active_set=[6,6,6,5],jlower=0,jhigher=2,exc=0,write_csf='rcsfmr.inp')
        zeros_6s5d= Rcsfgenerate('Ar',[AR_PREFIX + '4f(14,*)5d(1,*)6s(1,*)'],active_set=[6,5,5,4],jlower=0,jhigher=8,exc=0,write_csf='rcsfmr.inp')
        zeros_6p2= Rcsfgenerate('Ar',[AR_PREFIX + '4f(14,*)6p(2,*)'],active_set=[5,6,4,4],jlower=0,jhigher=6,exc=0,write_csf='rcsfmr.inp')
        zeros_even = zeros_6s2 + zeros_6s5d #+ zeros_6p2

        zeros_6s6p= Rcsfgenerate('Ar',[AR_PREFIX + '4f(14,*)6s(1,*)6p(1,*)',],active_set=[6,6,4,4],jlower=0,jhigher=4,exc=0,write_csf='rcsfmr.inp')
        zeros_4f5d= Rcsfgenerate('Ar',[AR_PREFIX + '4f(13,*)5d(1,*)6s(2,*)',],active_set=[6,5,5,4],jlower=0,jhigher=14,exc=0,write_csf='rcsfmr.inp')
        zeros_odd = zeros_6s6p + zeros_4f5d
        multireference = zeros_odd + zeros_even
        return multireference

    def ci_expansion_open_core(active_set,exc,write_csf,n_open = 6):
        # Only includes CV and VV!
        if n_open == 4:
            prefix = '3d(10,c)4s(2,*)4p(6,*)4d(10,*)5s(2,*)5p(6,*)'
        if n_open == 5:
            prefix = '3d(10,c)4s(2,i)4p(6,i)4d(10,i)5s(2,*)5p(6,*)'
        if n_open == 6:
            prefix = '3d(10,c)4s(2,i)4p(6,i)4d(10,i)5s(2,i)5p(6,i)'
        exp_6s2 = Rcsfgenerate('Ar',[prefix +  '4f(14,1)6s(2,1)'],active_set=active_set,jlower=0,jhigher=2,exc=exc,write_csf=write_csf)
        exp_4f5d = Rcsfgenerate('Ar',[prefix + '4f(14,1)5d(1,1)6s(1,1)'],active_set=active_set,jlower=0,jhigher=10,exc=exc,write_csf=write_csf)
        zeros_6p2= Rcsfgenerate('Ar',[prefix + '4f(14,i)6p(2,i)'],active_set=[5,6,4,4],jlower=exc,jhigher=6,exc=0,write_csf='rcsfmr.inp')
        exp_even = exp_6s2 + exp_4f5d #+ zeros_6p2

        exp_6s6p= Rcsfgenerate('Ar',[prefix + '4f(14,1)6s(1,1)6p(1,1)'],active_set=active_set,jlower=0,jhigher=6 ,exc=exc,write_csf=write_csf)
        exp_4f5d= Rcsfgenerate('Ar',[prefix + '4f(13,1)5d(1,1)6s(2,1)'],active_set=active_set,jlower=0,jhigher=16,exc=exc,write_csf=write_csf)
        exp_odd = exp_6s6p + exp_4f5d
        exp = exp_even + exp_odd
        return exp


    indices_mr = [[1],[1],[1],[1,2],[1],[1,2],[1],[1,2],[1],[1],[1],[1],[1],[1],[1]]
    def gen_expansion(calc_dir,active_set,exc,n_open):
        initialize(workdir=calc_dir,clist=clistordering)
        gen_multireference().execute(calc_dir)
        cmds = [
             Rnucleus(Z=70,A=172,neutralMass=171.936378,I=0,NDM=0,NQM=0),
             ci_expansion_open_core(active_set=active_set,exc = exc, write_csf = 'rcsf.inp',n_open = n_open),
             Rcsfinteract('Dirac-Coulomb'),
             Rangular(),
             Rwfnestimate(orbdict={'6p*':os.path.join(dir_6p, 'rwfn.out'),
                                   '5d*':os.path.join(dir_5d, 'rwfn.out'),
                                   '*' :os.path.join(dir_6s,'rwfn.out')},
                                   fallback='Thomas-Fermi'),
             Rmcdhf(indices_mr,orbs=[],specorbs=[],weighting_method='Standard',runs=1000),
             Rsave('safronova'),
             ]
        return [cmd.execute(workdir=calc_dir) for cmd in cmds]
    gen_expansion(dir_ci,active_set,exc,n_open)
    # input('Expansion OK?')
    Rci(calc_name='safronova',
         include_transverse=True,
         modify_freq=True,
         scale_factor='1.d-6',
         include_vacpol = True,
         include_nms= False,
         include_sms = False,
         est_self_energy= True,
         largest_n = 8,
         asfidx = indices_mr).execute(workdir = dir_ci)
    Rsave('safronova').execute(workdir = dir_ci)
    Rmixextract(calc_name = 'safronova',use_ci = True,tolerance = 0.05, sort = True, write_csf = 'rmix_5e-2.out').execute(dir_ci)
    Rlevels('safronova.cm').execute(dir_ci)
    #  run_ci(calc_dir = dir_cis[0],active_set = [7,7,6,5], exc = 1,n_open = 6)
    #  run_ci(calc_dir = dir_cis[1],active_set = [8,8,7,6], exc = 1,n_open = 6)
    #  run_ci(calc_dir = dir_cis[2],active_set = [9,9,8,7], exc = 1,n_open = 6)

    #  run_ci(calc_dir = dir_cis[3],active_set = [7,7,6,5], exc = 1,n_open = 5)
    #  run_ci(calc_dir = dir_cis[4],active_set = [8,8,7,6], exc = 1,n_open = 5)
    #  run_ci(calc_dir = dir_cis[5],active_set = [9,9,8,7], exc = 1,n_open = 5)

    #  run_ci(calc_dir = dir_cis[6],active_set = [7,7,6,5], exc = 1,n_open = 4)
    #  run_ci(calc_dir = dir_cis[7],active_set = [8,8,7,6], exc = 1,n_open = 4)
    #  run_ci(calc_dir = dir_cis[8],active_set = [9,9,8,7], exc = 1,n_open = 4)
import argparse
if __name__ == '__main__':
    parser = argparse.ArgumentParser('Neutral Yb MCDHF+CI script')
    parser.add_argument('--dir',help='Root directory to place calculation results', default = '/home/calvin/graspy/calc-outputs/safronova')
    #parser.add_argument('--maxn',help='Highest n of correlation orbital basis set',default=7,type=int)
    parser.add_argument('--maxn',help='Highest n of correlation orbital basis set', nargs='+', required=True,type= int)
    parser.add_argument('--coren',help='Lowest n of open core configurations',default = 6,type=int)
    parser.add_argument('--excitations',help='Number of excitations',default = 1,type=int)
    cmdargs = parser.parse_args()
    full_calculation(cmdargs.dir,active_set = cmdargs.maxn,exc = cmdargs.excitations,n_open = cmdargs.coren)
