from grasp import *
core_dir = '/home/calvin/graspy/calc-scripts/cores'
calc_dir = '../calc-outputs/ion-411-435'
initialize(workdir = calc_dir, clist =
['1s','2s','2p','3s','3p','3d','4s','4p','4d','5s','5p','4f','5d','6s','6p','5f','6d','7p','8s'])

#1. Yb core ([Xe]4f14) + 6s1 2S1/2, 5d1 2D3/2 and 2D5/2
# user-defined order of orbitals in clist.ref
# rnucleus_input="input_rnucleus_Z70"
# rnucleus < $rnucleus_input
# cat $rnucleus_input
# rcsfgenerate_input="input_6s_5d_rcsfgenerate"
# rcsfgenerate < $rcsfgenerate_input
# cat $rcsfgenerate_input
# cp rcsfgenerate.log 6s_5d.exc
# cp rcsf.out rcsf.inp
# printf 'y' | rangular
# rwfnestimate_input="input_6s_5d_rwfnestimate"
# rwfnestimate < $rwfnestimate_input
# cat $rwfnestimate_input
# rmcdhf_input="input_6s_5d_rmcdhf"
# rmcdhf < $rmcdhf_input
# cat $rmcdhf_input
# rsave 6s_5d
# save rwfn.out
# cp rwfn.out rwfn_6s_5d.out


def generate_6s_5d():
    cmdlist = [Rnucleus(Z=70,A=172,neutralMass=171.936378,I=0,NDM=0,NQM=0,rms_radius = 5.294,thickness = 2.18),
               Rcsfgenerate('Xe',['4f(14,c)6s(1,i)'],activeset=[6,5,4,4],jlower=1,jhigher=1,exc=0,ordering = 'User specified') + Rcsfgenerate('Xe',['4f(14,c)5d(1,i)'],activeset=[5,5,5,4],jlower=3,jhigher=5,exc=0,ordering = 'User specified'),
               Rangular(),
               Rwfnestimate(orbdict={'*' :os.path.join(core_dir,'rwfn_172Yb_core.out')}, fallback = 'Thomas-Fermi'),
               Rmcdhf(asfidx = [[1],[1],[1]],weightingmethod = 'Standard', orbs = ['*'], specorbs = ['*'],runs = 1000),
               Rsave('6s_5d')
               ]
    print(cmdlist)
    return [c.execute(workdir = calc_dir) for c in cmdlist]

generate_6s_5d()

def solve_orbs_2S12(active_levels,exc,active_set,optimize_orbs,rwfn_path,spec_orbs=[''],iccut = False,zero_first = False):
    cmdlist = []
    yb_inactive_core = '1s(2,c)2s(2,c)2p(6,c)3s(2,c)3p(6,c)3d(10,c)4s(2,c)4p(6,c)4d(10,c)'
    cmdlist.append(Rcsfgenerate('None',
                    csflist = [yb_inactive_core + active_levels],
                    activeset=active_set,
                    jlower=1,jhigher=1,exc=2))
    if zero_first:
        cmdlist += [Rcsfinteract('Dirac-Coulomb'),
                    Rcsfzerofirst(small_exp = f'../calc-scripts/cores/rcsf_2S12_mr0.99.inp',big_exp = 'rcsf.inp',write_csf = 'rcsf.inp')]
    if iccut is not False:
        cmdlist.append(Rangular(iccut = iccut))
    else:
        cmdlist.append(Rangular())
    cmdlist += [Rwfnestimate(orbdict={'*' : rwfn_path},
                             fallback = 'Thomas-Fermi'),
                Rmcdhf(asfidx = [[1]],
                       orbs = optimize_orbs,
                       specorbs=spec_orbs,runs = 100,
                       weightingmethod = 'Standard')]
    return [c.execute(workdir = calc_dir) for c in cmdlist]
solve_orbs_2S12(zero_first = False,
           active_levels = '5s(2,i)5p(6,i)4f(14,i)6s(1,i)',
           exc = 0,
           active_set = [6,5,4,4],
           optimize_orbs = ['*'],
           spec_orbs = ['*'],
           rwfn_path = '6s_5d.w') # take DHF rwfn
solve_orbs_2S12(zero_first = True,
               active_levels = '5s(2,*)5p(6,*)4f(14,*)6s(1,*)',
               exc = 2,
               active_set = [7,6,5,5,5],
               optimize_orbs = ['7s','6p-','6p','5d-','5d',
                          '5f-','5f','5g-','5g'],
               rwfn_path = os.path.join(calc_dir,'6s_5d.w')) # take DHF rwfn
solve_orbs_2S12(zero_first = True,
               active_levels = '5s(2,*)5p(6,*)4f(14,*)6s(1,*)',
               exc = 2,
               active_set = [8,7,6,6,5],
               optimize_orbs = ['8s','7p-','7p','6d-','6d',
                          '6f-','6f'],
               rwfn_path = os.path.join(calc_dir,'rwfn.out'))

Rsave('2S12').execute(workdir = calc_dir)
# Rci(calcname = '2S12',
#    includetransverse = True,
#    modifyfreq = True,
#    scalefactor = '1.d-6',
#    includevacpol = True,
#    includenms = False,includesms = False,estselfenergy = True,       largestn = 8,asfidx = [[1]])

# Rmixextract('2S12',useCI = False,tolerance = 0.01, sort = True).execute(workdir = calc_dir).execute(workdir)
# Rmixaccumulate('2S12',useCI=True,truncate_eps = 0.99, write_csf = 'rcsf.out').execute(workdir)

#### OLD SCRIPT ####
#1. separate 5d1 out from previous calculation (script_6s_5d)
#rnucleus_input="rnucleus_input_Z70"
#rnucleus < $rnucleus_input
#cat $rnucleus_input
#rcsfgenerate_input="input_2S12_rcsfgenerate_1"
#rcsfgenerate < $rcsfgenerate_input
#cp rcsf.out rcsf.inp
#cp rcsf.out rcsfmr.inp
#cat $rcsfgenerate_input
#TODO: Implement rcsfgenerate readout, e.g. cp rcsfgenerate.log 2S12.exc
#printf 'y' | rangular
#rwfnestimate_input="input_2S12_rwfnestimate_1"
#rwfnestimate < $rwfnestimate_input
#cat $rwfnestimate_input
#rmcdhf_input="input_2S12_rmcdhf_1"
#rmcdhf < $rmcdhf_input
#cat $rmcdhf_input


#2. add active orbital layer 1: 6s, 6p, 6d, 5f
#rcsfgenerate_input="input_2S12_rcsfgenerate_2"
#rcsfgenerate < $rcsfgenerate_input
#cat $rcsfgenerate_input
#cp rcsfgenerate.log 2S12.exc # TODO: READOUT
#cp rcsf.out rcsf.inp
#printf '1' | rcsfinteract
#cp rcsf.out rcsf.inp
## printf "rcsf_$stateID""_mr0.99.inp\nrcsf.inp" | rcsfzerofirst
#printf "rcsf_$stateID""_mr0.99.inp\nrcsf.inp" | rcsfzerofirst
#cp rcsf.out rcsf.inp
## printf 'y' | rangular
# printf 'y' | mpirun -np 4 rangular_mpi
#printf 'n\n147' | mpirun -np 4 rangular_mpi
#rwfnestimate_input="input_2S12_rwfnestimate_2"
#rwfnestimate < $rwfnestimate_input
#cat $rwfnestimate_input
#rmcdhf_input="input_2S12_rmcdhf_2"
# rmcdhf < $rmcdhf_input
#mpirun -np 4 rmcdhf_mpi < $rmcdhf_input
#cat $rmcdhf_input

##3. add active orbital layer 1: 7s, 7p, 7d, 6f
#rcsfgenerate_input="input_2S12_rcsfgenerate_3"
#rcsfgenerate < $rcsfgenerate_input
#cat $rcsfgenerate_input
#cp rcsfgenerate.log 2S12.exc
#cp rcsf.out rcsf.inp
#printf '1' | rcsfinteract
#cp rcsf.out rcsf.inp
##printf "rcsf_$stateID""_mr0.99.inp\nrcsf.inp" | rcsfzerofirst
##cp rcsf.out rcsf.inp
## printf 'y' | rangular#
#printf 'y' | mpirun -np 4 rangular_mpi
##printf 'n\n147' | mpirun -np 4 rangular_mpi
#rwfnestimate_input="input_2S12_rwfnestimate_3"
#rwfnestimate < $rwfnestimate_input
#cat $rwfnestimate_input
#rmcdhf_input="input_2S12_rmcdhf_3"
## rmcdhf < $rmcdhf_input
#mpirun -np 4 rmcdhf_mpi < $rmcdhf_input
#cat $rmcdhf_input


# rsave $stateID
# printf "1\n$stateID.w\n$stateID.w.readrwf" | readrwf
# printf "$stateID\nn\n0\ny" | rmixextract |& tee $stateID.rmix # DHF
# CI
#rci_input="input_2S12_rci"
#mpirun -np 4 rci_mpi < $rci_input
#cat $rci_input
# # jj2lsj
# jj2lsj_input="jj2lsj_input_2S12"
# jj2lsj < $jj2lsj_input
# cat $jj2lsj_input

# RIS4 & REDF
# printf "y\n$stateID\ny\ny\nn" | ris4
# printf "y\n$stateID\ny" | redf

# rmixextract & rmixaccumulate
# printf "$stateID\ny\n0.01\ny" | rmixextract |& tee $stateID.crmix # CI
# printf "$stateID\ny\n0.99\ny" | rmixaccumulate


# Rlevels(['2S12.m','2Dduplet.m'])
# rlevelseV $stateID.m
# rlevelseV $stateID.m 2Dduplet.m
# rlevelseV $stateID.cm 2Dduplet.cm
