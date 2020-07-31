from graspy.grasp import *
import sys
from shutil import copyfile
from os import path

testdir = './calc-outputs/radium2out'
initialize(workdir=testdir)

def gen_nucleus(calc_dir):
    cmds = [Rnucleus(Z=88,A=226,neutralMass=225.95172496853601 ,I=0,NDM=0,NQM=0,rms_radius = 5.6622426556969394,thickness = 2.2999999999999998)]
    return [cmd.execute(workdir = calc_dir) for cmd in cmds]

def angular_integration(calc_dir):
    return Rangular().execute(workdir = calc_dir)

def estimate_wavefunctions(calc_dir,previous_rwfn = None):
    if previous_rwfn is None:
        return Rwfnestimate(orbdict = None,fallback =  'Thomas-Fermi').execute(workdir=testdir)
    else:
        return Rwfnestimate(orbdict = {'*': previous_rwfn},fallback = 'Thomas-Fermi').execute(workdir = calc_dir)
    

def run_hartree_fock_save(calc_dir,grid = None,orbs = ['*'],specorbs = ['*'],integration_method = None,calc_name = None):
    out = Rmcdhf(asfidx = [[1],[1,2],[1,2],[1,2,3,4],[1,2,3,4],[1,2],[1,2]],orbs = orbs, specorbs = specorbs, runs = 20000, weighting_method = 'Standard',integration_method = integration_method).execute(workdir = calc_dir)
    if calc_name is not None:
        Rsave(calc_name).execute(workdir = calc_dir)
        return os.path.join(calc_dir,calc_name) + '.w'
    else:
        return out

gen_nucleus(testdir)
ref_7s = Rcsfgenerate(core='Rn',ordering = 'Default',
            csflist=['7s(2,i)'],
            active_set=[7,6,5,4],
            jlower=0,jhigher=0,exc=0, write_csf= 'rcsfmr.inp')
ref_7p = Rcsfgenerate(core='Rn',ordering = 'Default',
            csflist=['7s(1,i)7p(1,i)'],
            active_set = [7,7,5,4],
            jlower=0,jhigher=4,exc=0,write_csf= 'rcsfmr.inp')
ref6_d = Rcsfgenerate(core='Rn',ordering = 'Default',
            csflist=['7s(1,i)6d(1,i)'],
            active_set = [7,6,6,4],
            jlower=2,jhigher=6,exc=0,write_csf= 'rcsfmr.inp')
            
mr = ref_7p + ref_7s + ref6_d
mr.execute(workdir=testdir)

angular_integration(testdir)

principle_quantum = [1, 2, 3, 4, 5, 6 ,7]

orb_list = [['1s'],         #list of orbitals to be passed into each n iteration
['1s','2s*','2p*'],
['1s','2s*','2p*','3s','3p*','3d*'],
['1s','2s*','2p*','3s','3p*','3d*','4s','4p*','4d*', '4f*'],
['1s','2s*','2p*','3s','3p*','3d*','4s','4p*','4d*', '4f*', '5s','5p*','5d*','5f*'],
['1s','2s*','2p*','3s','3p*','3d*','4s','4p*','4d*', '4f*', '5s','5p*','5d*','5f*', '6s', '6p*','6d*'],
['1s','2s*','2p*','3s','3p*','3d*','4s','4p*','4d*', '4f*', '5s','5p*','5d*','5f*', '6s', '6p*','6d*','7s','7p*']]

old_rwfn = None
for orb_config in orb_list[0:1]:
    print("Running ", *(orb_config))
    estimate_wavefunctions(testdir,old_rwfn)
    old_rwfn = run_hartree_fock_save(testdir, grid=None,orbs = orb_config,specorbs = ['*'],integration_method = None,calc_name = None)
    
    
    
    #rwfn_cp_file = 'rwfn' + str(i) + '.inp'
    #copyfile(path.join(testdir, 'rwfn.inp'),path.join(testdir, rwfn_cp_file))
    #old_rwfn=rwfn_cp_file
#run_hartree_fock_save(testdir, grid=None, orbs=['*'],specorbs = ['*'],integration_method = None,calc_name = None)

