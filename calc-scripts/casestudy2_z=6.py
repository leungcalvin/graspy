from grasp import *
import shutil

testdir = '/Users/admin/Desktop/asdrp2020/grasp/grasptest/case2/script/Z6'
initialize(workdir=testdir)

nucleus = Rnucleus(Z = 6, A = 12, neutralMass = 12.0107, I = 1, NDM = 1, NQM = 1)
nucleus.execute(workdir=testdir)

ref_2s = Rcsfgenerate(core='None',ordering = 'Default',
            csflist=['1s(2,i)2s(1,i)'],
            active_set=[2,2],
            jlower=1,jhigher=1,exc=0, write_csf= 'DF.c')
ref_2p = Rcsfgenerate(core='None',ordering = 'Default',
            csflist=['1s(2,i)2p(1,i)'],
            active_set=[2,2],
            jlower=1,jhigher=3,exc=0, write_csf= 'DF.c')

mr = ref_2s + ref_2p
mr.execute(workdir = testdir)

# 2. Generate expansions
even_exp = Rcsfgenerate(core='None',ordering = 'Default',
            csflist=['1s(2,*)2s(1,*)'],
            active_set=[5,5,5,5,5],
            jlower=1,jhigher=1,exc=3, write_csf='even.c')

even_exp.execute(workdir=testdir)

even_exp_split = Rcsfsplit(calc_name = 'even',nsplit = 3, splitorbs = range(3,5))
even_exp_split.execute(workdir = testdir)

odd_exp = Rcsfgenerate(core='None',ordering = 'Default',
            csflist=['1s(2,*)2p(1,*)'],
            active_set=[5,5,5,5,5],
            jlower=1,jhigher=3,exc=3, write_csf= 'odd.c')

odd_exp.execute(workdir=testdir)

odd_exp_split = Rcsfsplit(calc_name = 'odd',nsplit=3,splitorbs = range(3,5))
odd_exp_split.execute(workdir = testdir)

# 3. Ground and excited reference states
# for z in range(6,12):
    # cp ../DF.c rcsf.inp
path_to_copy = os.path.join(testdir,'DF.c')
target_path = os.path.join(testdir,'rcsf.inp')
shutil.copyfile(path_to_copy,target_path)
ref_states_1s2s = [Rangular(), 
    # Get angular data
    Rwfnestimate(orbdict = None, fallback='Thomas-Fermi'), 
    #Get initial estimates of wave functions
    Rmcdhf([[1],[1],[1]],orbs = ['*'],specorbs = ['*'], runs = 100, weighting_method = 'Standard'), 
    # Perform self-consistent field calculations
    #output =  outodd_rmcdhf_initial 
    Rsave('DF')]
[cmd.execute(workdir = testdir) for cmd in ref_states_1s2s]
# cp DF.w even2.w
path_to_copy_1 = os.path.join(testdir,'DF.w')
target_path_1 = os.path.join(testdir,'even2.w')
shutil.copyfile(path_to_copy_1,target_path_1)
# cp DF.w odd2.w
path_to_copy_2 = os.path.join(testdir,'DF.w')
target_path_2 = os.path.join(testdir,'odd2.w')
shutil.copyfile(path_to_copy_2,target_path_2)


# 4. Perform calculations for the even states
for n in range(3,6):
    # for z in range(6,12):
    # cp ../even${n}.c rcsf.inp
    path_to_copy = os.path.join(testdir,f'even{n}.c')
    target_path = os.path.join(testdir,'rcsf.inp')
    shutil.copyfile(path_to_copy,target_path)
    k = n - 1
    print(k)
    calculations_even = [Rangular(), 
        # Get angular data
        Rwfnestimate(orbdict = {'*':f'even{k}.w'},fallback = 'Thomas-Fermi'), #how to represent multiple files with differing numbers? 
        # Get initial estimates of wave functions
        Rmcdhf(asfidx = [[1]],orbs = ['*'],specorbs = ['*'], runs = 100, weighting_method = 'Standard'), 
        # Perform self-consistent field calculations
        #Output = outodd_rmcdhf_${n}
        Rsave(f'even{n}')]
    [cmd.execute(workdir = testdir) for cmd in calculations_even]
    n += 1

# 5. Perform calculations for the odd states
for n in range(3,6):
    # for z in range(6,12):
    # cp ../odd${n}.c rcsf.inp
    path_to_copy_3 = os.path.join(testdir,f'odd{n}.c')
    target_path_3 = os.path.join(testdir,'rcsf.inp')
    shutil.copyfile(path_to_copy_3,target_path_3)
    k = n - 1
    print(k)
    calculations_odd = [Rangular(), 
            # Get angular data
            Rwfnestimate(orbdict = {'*':f'odd{k}.w'},fallback = 'Thomas-Fermi'), #how to represent multiple files with differing numbers? 
            # Get initial estimates of wave functions
            Rmcdhf(asfidx = [[1],[1]],orbs = ['*'],specorbs = ['*'], runs = 100, weighting_method = 'Standard'), 
            # Perform self-consistent field calculations
            # output = outodd_rmcdhf_${n}
            Rsave(f'odd{n}')]
    [cmd.execute(workdir = testdir) for cmd in calculations_odd]
    n += 1

# 6. Configuration interaction and transition calculations
# for z in range(6,12):
calculations_even5 = Rci(calc_name= 'even5',
    include_transverse=True,
    modify_freq=True,
    scale_factor='1.d-6',
    include_vacpol = True,
    include_nms= False,
    include_sms = False,
    est_self_energy= True,
    largest_n = 3,
    asfidx = [[1]])
calculations_even5.execute(workdir=testdir)
# output = outeven_rci 
calculations_odd5 = Rci(calc_name= 'odd5',
    include_transverse=True,
    modify_freq=True,
    scale_factor='1.d-6',
    include_vacpol = True,
    include_nms= False,
    include_sms = False,
    est_self_energy= True,
    largest_n = 3,
    asfidx = [[1],[1]])
calculations_odd5.execute(workdir=testdir)

save_Ang_data = [Rbiotransform(use_ci=True,calc_name_initial = 'even5',calc_name_final = 'odd5', transform_all = True), 
# Run rbiotransform ans save angular data
    Rtransition(use_ci=True,calc_name_initial = 'even5',calc_name_final = 'odd5',transition_spec = ['E1'])]
# Run rtransition save angular data        
[cmd.execute(workdir = testdir) for cmd in save_Ang_data]

