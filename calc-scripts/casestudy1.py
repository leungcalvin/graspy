from grasp import *
import shutil
#
# This implements the GRASP 2018 script demonstration as found in the GRASP 2018 manual, using the GRASPy interface.
#
testdir = '/Users/admin/Desktop/asdrp2020/grasp/grasptest/case1/script'
initialize(workdir=testdir)

#Generate odd CSF expansions  1.1 MR for 2s(2)2p, 2p(3)
odd_mr = Rcsfgenerate(core='None',ordering = 'Default',
            csflist=['1s(2,i)2s(2,i)2p(1,i)',                                                                                                                                                                                                                                                                                       
                     '1s(2,i)2p(3,i)'],
            active_set=[2,2],
            jlower=1,jhigher=3,exc=0, write_csf= 'odd2.c') 
            #cp rcsf.out odd2.c -> do i need to add this? does it copy over?

odd_exp = Rcsfgenerate(core='None',ordering = 'Default',
            csflist=['1s(2,*)2s(2,*)2p(1,*)',
                     '1s(2,*)2p(3,*)'],
            active_set=[6,6,6,6,6,6],
            jlower=1,jhigher=3,exc=2, write_csf= 'odd.c')
            #cp rcsf.out odd.c -> do i need to add this? does it copy over?

odd_mr.execute(workdir = testdir) #the execute() method performs the actual call to GRASP 2018.

# # Split into odd3.c, odd4.c, odd5.c, odd6.c
# odd_exp_split = Rcsfsplit(calc_name = 'odd',nsplit=4,splitorbs = range(3,7))
# odd_exp_split.execute(workdir = testdir)

# #Generate even CSF expansions 2.1 for 2s2p(2)
# even_2s2p2_mr = Rcsfgenerate(core='None',ordering = 'Default',
#             csflist=['1s(2,i)2s(1,i)2p(2,i)',
#                      ],
#             active_set=[2,2],
#             jlower=1,jhigher=5,exc=0, write_csf= 'even2.c')
#             #cp rcsf.out even2.c -> do i need to add rcsf.out?

# even_2s2p2_exp = Rcsfgenerate(core='None',ordering = 'Default',
#             csflist=['1s(2,*)2s(1,*)2p(2,*)',
#                      ],
#             active_set=[6,6,6,6,6,6],
#             jlower=1,jhigher=5,exc=2, write_csf='even.c')
#             #cp rcsf.out even.c -> do i need to add rcsf.out?

# even_2s2p2_mr.execute(workdir = testdir)

# #Split into even3.c, even4.c, even5.c, even6.c
# even_2s2p2_exp_split = Rcsfsplit(calc_name = 'even',nsplit = 4, splitorbs = range(3,7))
# even_2s2p2_exp_split.execute(workdir = testdir)

# #2. Get Nuclear Data
# nuclear_data = Rnucleus(Z = 42, A = 96, neutralMass = 96, I = 1, NDM = 1,NQM = 1)
# nuclear_data.execute(workdir = testdir)

# #3. Get initial estimate
# #For n=2, Get initial estimates for odd.
# #copy odd2.c to rcsf.inp
# path_to_copy = os.path.join(testdir,'odd2.c')
# target_path = os.path.join(testdir,'rcsf.inp')
# shutil.copyfile(path_to_copy,target_path)

# init_est_odd = [Rangular(), #Get initial estimates of wave functions
#         Rwfnestimate(orbdict = None, fallback='Thomas-Fermi'), 
#         Rmcdhf([[1],[1],[5]],orbs = ['*'],specorbs = ['*'], runs = 100, weighting_method = 'Standard'), 
#         # Perform self-consistent field calculations
#         #output =  outodd_rmcdhf_initial 
#         Rsave('odd2.c')]
# [cmd.execute(workdir = testdir) for cmd in init_est_odd]

# #For n=2, Get initial estimates for even
# #copy even2.c to rcsf.inp
# path_to_copy = os.path.join(testdir,'even2.c')
# target_path = os.path.join(testdir,'rcsf.inp')
# shutil.copyfile(path_to_copy,target_path)

# inti_est_even = [Rangular(), 
#         Rwfnestimate(orbdict = None, fallback='Thomas-Fermi'),
#         # Get initial estimates of wave functions
#         Rmcdhf([[1,2,3],[1,2,3],[1,2], [5]],orbs = ['*'],specorbs = ['*'], runs = 100, weighting_method = 'Standard'), 
#         # Perform self-consistent field calculations
#         # output = outeven_rmcdhf_initial
#         Rsave('even2.c')]
# [cmd.execute(workdir = testdir) for cmd in inti_est_even]

#4. rmcdhf and rci calculation
#Get results for odd n=3,4,5,6
for n in range(3,6):
        path_to_copy_2 = os.path.join(testdir,f'odd{n}.c') 
        target_path_2 = os.path.join(testdir,'rcsf.inp')
        shutil.copyfile(path_to_copy_2,target_path_2)
        m = n - 1
        print(m)
        print(n)
        calculations_odd = [Rangular(), 
                # Get angular data
                Rwfnestimate(orbdict = {'*':f'odd{m}.w'},fallback = 'Thomas-Fermi'), #how to represent multiple files with differing numbers?
                # Get initial estimates of wave functions
                Rmcdhf(asfidx = [[1],[1],[5]],orbs = ['*'],specorbs = ['*'], runs = 100, weighting_method = 'Standard'), 
                #output = outodd_rmcdhf_${n}
                # Perform self-consistent field calculations
                Rsave(f'odd{n}')] #rsave odd${n}
                #should it be f'odd{n}' or 'odd'
        [cmd.execute(workdir = testdir) for cmd in calculations_odd]
        n += 1

        ## Perform Breit-correction using RCI for n=6. First copy to other file names
        if n == 6:
                #cp odd${n}.c oddCI${n}.c
                path_to_copy_3 = os.path.join(testdir,'odd6.c')
                target_path_3 = os.path.join(testdir,'oddCI6.c')
                shutil.copyfile(path_to_copy_3,target_path_3)
                #cp odd${n}.w oddCI${n}.w
                path_to_copy_4 = os.path.join(testdir,'odd6.w')
                target_path_4 = os.path.join(testdir,'oddCI6.w')
                shutil.copyfile(path_to_copy_4,target_path_4)
                BC_odd = Rci(calc_name= 'oddCI6', #oddCI${n}
                #should it be f'oddCI{n}' or 'oddCI6'
                        include_transverse=True,
                        modify_freq=True,
                        scale_factor='1.d-6',
                        include_vacpol = True,
                        include_nms= False,
                        include_sms = False,
                        est_self_energy= True,
                        largest_n = 4,
                        asfidx = [[1],[1]]) 
                        #output = outodd_rci
                BC_odd.execute(workdir=testdir)
        else:
                continue
        # transform to LSJ-coupling
        Transform_LSJ_odd = JJtoLSJ(calc_name= f'oddCI{n}',use_ci = True, unique = True)
        #should it be f'oddCI{n}' or 'oddCI'
        Transform_LSJ_odd.execute(workdir = testdir)

# Get results for even n=3,4,5,6
for n in range(3,6):
        #cp even${n}.c rcsf.inp
        path_to_copy_5 = os.path.join(testdir,f'even{n}.c') #how to represent multiple files with differing numbers?
        target_path_5 = os.path.join(testdir,'rcsf.inp')
        shutil.copyfile(path_to_copy_5,target_path_5)
        m = n -1
        print(m)
        print(n)
        calculations_even = [Rangular(), 
                # Get angular data
                Rwfnestimate(orbdict = {'*':f'even{m}.w'},fallback = 'Thomas-Fermi'), #how to represent multiple files with differing numbers? 
                # Get initial estimates of wave functions
                Rmcdhf(asfidx = [[1,2,3],[1,2,3],[1,2],[5]],orbs = ['*'],specorbs = ['*'], runs = 100, weighting_method = 'Standard'), 
                # Perform self-consistent field calculations
                #Output = outodd_rmcdhf_${n}
                Rsave(f'even{n}')]
                #should it be f'even{n}' or 'even'
        [cmd.execute(workdir = testdir) for cmd in calculations_even]
        n += 1

        # Perform Breit-correction using RCI for n=6
        if n == 6:
                #cp even${n}.c evenCI${n}.c
                path_to_copy_6 = os.path.join(testdir,'even6.c') 
                target_path_6 = os.path.join(testdir,'evenCI6.c') 
                shutil.copyfile(path_to_copy_6,target_path_6)
                #cp even${n}.w evenCI${n}.w
                path_to_copy_7 = os.path.join(testdir,'even6.w') 
                target_path_7 = os.path.join(testdir,'evenCI6.w') 
                shutil.copyfile(path_to_copy_7,target_path_7)
                BC_even = Rci(calc_name= 'evenCI6',
                #should it be f'evenCI{n}' or 'evenCI6'
                        include_transverse=True,
                        modify_freq=True,
                        scale_factor='1.d-6',
                        include_vacpol = True,
                        include_nms= False,
                        include_sms = False,
                        est_self_energy= True,
                        largest_n = 4,
                        asfidx = [[1,2,3],[1,2,3],[1,2],[5]])
                        #output = outeven_rci
                BC_even.execute(workdir=testdir)
        else: 
                continue
        # transform to LSJ-coupling
        Transform_LSJ_even = JJtoLSJ(calc_name= f'evenCI{n}',use_ci = True, unique = True)
        #should it be f'evenCI{n}' or 'evenCI'
        Transform_LSJ_even.execute(workdir = testdir) 

# 5. Transition Calculation
# Perform transition calculation for the n=6 CI results
# n = 6
BT = [Rbiotransform(use_ci=True,calc_name_initial = 'oddCI6',calc_name_final = 'evenCI6', transform_all = True), 
        #should it be f'oddCI{n}' or 'oddCI', #should it be f'evenCI{n}' or 'evenCI'
        #First the biorthonormal transformations
        Rtransition(use_ci=True,calc_name_initial = 'oddCI6',calc_name_final = 'evenCI6',transition_spec = ['E1'])]
        #should it be f'oddCI{n}' or 'oddCI', #should it be f'evenCI{n}' or 'evenCI'
        #Then the transition calculations
[cmd.execute(workdir = testdir) for cmd in BT]
