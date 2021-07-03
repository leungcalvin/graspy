from graspy.grasp import *
import sys
from shutil import copyfile
from os import path

# default n (open core orbitals) and excitations.
princ_q_num = 10
excitations = 1

if len(sys.argv) > 1 and sys.argv[1] != None :
    princ_q_num = int(sys.argv[1])

if len(sys.argv) > 2 and sys.argv[2] != None :
    excitations = int(sys.argv[2])

print("Running for n =", princ_q_num )
print("Running with Exications =", excitations )



#testdir = './calc-outputs/radium2out'
testdir = './calc-outputs/radium_' + str(princ_q_num) + '_' + str(excitations)
print("Using testdir =", testdir)



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
    out = Rmcdhf(asfidx = [[1],[1],[1],[1,2],[1,2],[1,2,3,4,5],[1],[1,2,3],[1]],orbs = orbs, specorbs = specorbs, runs = 20000, weighting_method = 'Standard',integration_method = integration_method).execute(workdir = calc_dir)
    Rsave(calc_name).execute(workdir = calc_dir)
    return calc_name + '.w'
    #else:
        #return out

#[1],[1],[1],[1,2],[1,2],[1],[1]
#[1],[1],[1,2],[1,2],[1,2,3,4,5],[1],[1,2,3],[1]
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
ref_7p_6d = Rcsfgenerate(core='Rn',ordering = 'Default',
            csflist=['7p(1,i)6d(1,i)'],
            active_set = [6,7,6,4],
            jlower=4,jhigher=8,exc=0,write_csf= 'rcsfmr.inp')
            
mr = ref_7p + ref_7s + ref6_d + ref_7p_6d
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

# set the right active sets based on n
if princ_q_num == 8:
    a_set = [8,8,8,8,8]
elif princ_q_num == 10:
    a_set = [10,10,10,10,10]
elif princ_q_num == 11:
    a_set = [11,11,11,11,11,11]
elif princ_q_num == 12:
    a_set = [12,12,12,12,12,12]

ref_7s_e = Rcsfgenerate(core='Xe',ordering = 'Default',
            csflist=['4f(14,i)5d(10,9)6s(2,1)6p(6,5)7s(2,*)'],
            #active_set=[10,10,10,10,10],
            active_set=a_set,
            jlower=0,jhigher=0,exc=excitations, write_csf= 'rcsf.inp')
ref_7p_e = Rcsfgenerate(core='Xe',ordering = 'Default',
            csflist=['4f(14,i)5d(10,9)6s(2,1)6p(6,5)7s(1,*)7p(1,*)'],
            #active_set = [10,10,10,10,10],
            active_set=a_set,
            jlower=0,jhigher=4,exc=excitations,write_csf= 'rcsf.inp')
ref6_d_e = Rcsfgenerate(core='Xe',ordering = 'Default',
            csflist=['4f(14,i)5d(10,9)6s(2,1)6p(6,5)7s(1,*)6d(1,*)'],
            #active_set = [10,10,10,10,10],
            active_set=a_set,
            jlower=2,jhigher=6,exc=excitations,write_csf= 'rcsf.inp')
ref_7p_6d_e = Rcsfgenerate(core='Xe',ordering = 'Default',
            csflist=['4f(14,i)5d(10,9)6s(2,1)6p(6,5)7s(2,1)7p(1,*)6d(1,*)'],
            #active_set = [10,10,10,10,10],
            active_set=a_set,
            jlower=4,jhigher=8,exc=excitations,write_csf= 'rcsf.inp')

mr_e = ref_7p_e + ref_7s_e + ref6_d_e + ref_7p_6d_e
mr_e.execute(workdir=testdir)

angular_integration(testdir)
estimate_wavefunctions(testdir, old_rwfn, fallback='Screened Hydrogenic')
old_rwfn = run_hartree_fock_save(testdir, grid=None,orbs=[], specorbs = [],integration_method = 3,calc_name = 'eighth')
print("old_rwfn=",old_rwfn)
# saves eight.w and eight.c into radium.w and radium.c respectively
Rsave('radium').execute(workdir = testdir)

# RCI uses radium.c and radium.w and produces radium.cm, radium.clog, radium.m,rci.res,  radium.csum...
Rci(calc_name= 'radium',
    include_transverse=True,
    modify_freq=True,
    scale_factor='1.d-6',
    include_vacpol = True,
    include_nms= False,
    include_sms = False,
    est_self_energy= True,
    largest_n = 8,
    asfidx = [[1],[1],[1],[1,2],[1,2],[1,2,3,4,5],[1],[1,2,3],[1]]).execute(workdir = testdir)

#Rsave('radiumCI').execute(workdir = testdir)
# rlevels uses radium.cm to work. It is not happy with just radium
Rlevels('radium.cm').execute(workdir = testdir)

