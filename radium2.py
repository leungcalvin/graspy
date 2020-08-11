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

def estimate_wavefunctions(calc_dir,previous_rwfn = None,fallback='Thomas-Fermi'):
    if previous_rwfn is None:
        return Rwfnestimate(orbdict = None,fallback = fallback).execute(workdir=testdir)
    else:
        return Rwfnestimate(orbdict = {'*': previous_rwfn},fallback =fallback).execute(workdir = calc_dir)
    

def run_hartree_fock_save(calc_dir,grid = None,orbs = ['*'],specorbs = ['*'],integration_method = None,calc_name = None):
    out = Rmcdhf(asfidx = [[1],[1],[1],[1,2],[1,2],[1],[1]],orbs = orbs, specorbs = specorbs, runs = 20000, weighting_method = 'Standard',integration_method = integration_method).execute(workdir = calc_dir)
    Rsave(calc_name).execute(workdir = calc_dir)
    return calc_name + '.w'
    #else:
        #return out

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


orb_list = [['1*'], 
['1*', '2*'],
['1*', '2*', '3*'], 
['1*', '2*', '3*', '4*'],
['1*', '2*', '3*', '4*', '5*'],
['1*', '2*', '3*', '4*','5*','6s*','6p*','6d*'],
['1*', '2*', '3*', '4*','5*','6s*','6p*','6d*','7s*','7p*']]

'''

orb_config_dict = {'first':['1*'],      #dictionary with key=calcdir and value=orbital because these are the two variables in the loop
'second':['1*', '2*'], 
'third':['1*', '2*', '3*'], 
'fourth':['1*', '2*', '3*', '4*'], 
'fifth':['1*', '2*', '3*', '4*', '5*'], 
'sixth':['1*', '2*', '3*', '4*','5*','6s*','6p*','6d*'],
'seventh':['1*', '2*', '3*', '4*','5*','6s*','6p*','6d*','7s*','7p*']}


old_rwfn = None     #created a loop that takes the corresponding dictionary items and keeps revising each configuration
for config_item in orb_config_dict:
    print("key=",config_item,"value=",orb_config_dict[config_item])
    estimate_wavefunctions(testdir, old_rwfn)
    print("old_rwfn=",old_rwfn)
    old_rwfn = run_hartree_fock_save(testdir, grid=None,orbs=orb_config_dict[config_item], specorbs = ['*'],integration_method = 3,calc_name = config_item)
'''

#ran each orbital manually 
 
old_rwfn = None
estimate_wavefunctions(testdir,old_rwfn)
old_rwfn = run_hartree_fock_save(testdir, grid=None,orbs= orb_list[0], specorbs = ['*'],integration_method = 3,calc_name = 'first')
print("old_rwfn=",old_rwfn)
print("testdir = ",testdir)

estimate_wavefunctions(testdir,old_rwfn)
old_rwfn = run_hartree_fock_save(testdir, grid=None,orbs=orb_list[1], specorbs = ['*'],integration_method = 3,calc_name = 'second')
print("old_rwfn=",old_rwfn)

estimate_wavefunctions(testdir,old_rwfn)
old_rwfn = run_hartree_fock_save(testdir, grid=None,orbs=orb_list[2], specorbs = ['*'],integration_method = 3,calc_name = 'third')
print("old_rwfn=",old_rwfn)

estimate_wavefunctions(testdir,old_rwfn)
old_rwfn = run_hartree_fock_save(testdir, grid=None,orbs=orb_list[3], specorbs = ['*'],integration_method = 3,calc_name = 'fourth')
print("old_rwfn=",old_rwfn)

estimate_wavefunctions(testdir,old_rwfn)
old_rwfn = run_hartree_fock_save(testdir, grid=None,orbs=orb_list[4], specorbs = ['*'],integration_method = 3,calc_name = 'fifth')
print("old_rwfn=",old_rwfn)

estimate_wavefunctions(testdir,old_rwfn)
old_rwfn = run_hartree_fock_save(testdir, grid=None,orbs=orb_list[5], specorbs = ['*'],integration_method = 3,calc_name = 'sixth')
print("old_rwfn=",old_rwfn)

estimate_wavefunctions(testdir,old_rwfn)
old_rwfn = run_hartree_fock_save(testdir, grid=None,orbs=orb_list[6], specorbs = ['*'],integration_method = 3,calc_name = 'seventh')
print("old_rwfn=",old_rwfn)


ref_7s_e = Rcsfgenerate(core='Rn',ordering = 'Default',
            csflist=['7s(2,*)'],
            active_set=[8,8,8,8,8],
            jlower=0,jhigher=0,exc=1, write_csf= 'rcsf.inp')
ref_7p_e = Rcsfgenerate(core='Rn',ordering = 'Default',
            csflist=['7s(1,*)7p(1,*)'],
            active_set = [8,8,8,8,8],
            jlower=0,jhigher=4,exc=1,write_csf= 'rcsf.inp')
ref6_d_e = Rcsfgenerate(core='Rn',ordering = 'Default',
            csflist=['7s(1,*)6d(1,*)'],
            active_set = [8,8,8,8,8],
            jlower=2,jhigher=6,exc=1,write_csf= 'rcsf.inp')
            
mr_e = ref_7p_e + ref_7s_e + ref6_d_e
mr_e.execute(workdir=testdir)

angular_integration(testdir)
estimate_wavefunctions(testdir, old_rwfn, fallback='Screened Hydrogenic')
old_rwfn = run_hartree_fock_save(testdir, grid=None,orbs=[], specorbs = [],integration_method = 3,calc_name = 'eighth')
print("old_rwfn=",old_rwfn)
Rsave('radium').execute(workdir = testdir)

Rci(calc_name= 'radium',
    include_transverse=True,
    modify_freq=True,
    scale_factor='1.d-6',
    include_vacpol = True,
    include_nms= False,
    include_sms = False,
    est_self_energy= True,
    largest_n = 8,
    asfidx = [[1],[1],[1],[1,2],[1,2],[1],[1]]).execute(workdir = testdir)

#Rsave('radiumCI').execute(workdir = testdir)
#Rlevels('radiumCI').execute(workdir = testdir)