from graspy.grasp import *
#
# This implements the GRASP 2018 calculation for 1s2 2s 2S and 1s2 2p 2P in Li I, as found in the GRASP 2018 manual, using the GRASPy interface.
#
ex1_workdir = '/home/cleu/example1_mpi'
initialize(workdir=ex1_workdir)

# 1) Generate a multireference consisting of the 1s2 2s and 1s2 2p configurations
ref_2s = Rcsfgenerate(core='None',ordering = 'Default',
            csflist=['1s(2,i)2s(1,i)'],
            active_set=[2],
            jlower=1,jhigher=1,exc=0,write_csf= 'rcsfmr.inp')
ref_2p = Rcsfgenerate(core='None',ordering = 'Default',
            csflist=['1s(2,i)2p(1,i)'],
            active_set = [1,2],
            jlower=1,jhigher=3,exc=0,write_csf= 'rcsfmr.inp')

# The three states which make up the multireference can be added together with arithmetic.
mr = ref_2s + ref_2p
mr.execute(workdir = ex1_workdir) #the execute() method performs the actual call to GRASP 2018.
# 2) Perform a SCF procedure to solve for the 1s,2s,2p orbitals.
MR_DHF =[
        Rnucleus(Z=3,A=7,neutralMass=6.941,I=1.5,NDM=3.2564268,NQM=-0.04),
        Rangular(),
        Rwfnestimate(orbdict = None, fallback='Thomas-Fermi'),
        ]
for cmd in MR_DHF:
    cmd.execute(workdir = ex1_workdir)

Rmcdhf([[1],[1],[1]],orbs = ['*'],specorbs = ['*'], runs = 100, weighting_method = 'Standard').execute(workdir = ex1_workdir)
Rsave('2s_2p_DF').execute(workdir = ex1_workdir)

# 3) Generate a CAS expansion from the 2S configuration.
CAS_2S_exp = Rcsfgenerate(core='None',ordering = 'Default',
            csflist=['1s(2,*)2s(1,*)'],
            active_set=[3,3,3],
            jlower=1,jhigher=1,exc=3,write_csf = 'rcsf.inp')
CAS_2S_exp.execute(workdir = ex1_workdir)

# 4) Solve for the n=3 correlation orbitals, using orbitals generated from 2s_2p_DF.w above.
indices_2S = [[1]]
Rangular().execute(workdir = ex1_workdir)
Rwfnestimate(orbdict = {'*':'2s_2p_DF.w'},fallback = 'Thomas-Fermi').execute(workdir = ex1_workdir)
Rmcdhf(asfidx = indices_2S,
    orbs = ['3*'],
    specorbs = [' '],
    runs = 100, weighting_method = 'Standard').execute(workdir = ex1_workdir)
Rsave('2s_3').execute(workdir = ex1_workdir)

# 5) Perform CI on the 2S expansion.
Rci(calc_name='2s_3',
    include_transverse=True,
    modify_freq=True,
    scale_factor='1.d-6',
    include_vacpol = True,
    include_nms= False,
    include_sms = False,
    est_self_energy= True,
    largest_n = 3,
    asfidx = indices_2S).execute(workdir = ex1_workdir),

JJtoLSJ(calc_name= '2s_3',use_ci = True, unique = True).execute(workdir = ex1_workdir)

# 6) Generate a CAS expansion from the 2P configuration.
CAS_2P_exp = Rcsfgenerate(core='None',ordering = 'Default',
            csflist=['1s(2,*)2p(1,*)'],
            active_set=[3,3,3],
            jlower=1,jhigher=3,exc=3)

CAS_2P_exp.execute(workdir = ex1_workdir)

# 7) Solve for the n=3 correlation orbitals, using orbitals generated from 2s_2p_DF.w above.
indices_2P = [[1],[1]]
Rangular().execute(workdir = ex1_workdir)
Rwfnestimate(orbdict = {'*':'2s_2p_DF.w'},fallback = 'Thomas-Fermi').execute(workdir = ex1_workdir)
Rmcdhf(asfidx = indices_2P,orbs = ['3*'],specorbs = [' '],runs = 100, weighting_method = 'Standard').execute(workdir = ex1_workdir)
Rsave('2p_3').execute(workdir = ex1_workdir)

# 8) Perform CI on the 2P expansion.

Rci(calc_name='2p_3',
    include_transverse=True,
    modify_freq=True,
    scale_factor='1.d-6',
    include_vacpol = True,
    include_nms= False,
    include_sms = False,
    est_self_energy= True,
    largest_n = 3,
    asfidx = indices_2P).execute(workdir = ex1_workdir),
JJtoLSJ(calc_name= '2p_3',use_ci = True, unique = True).execute(workdir = ex1_workdir)

# 9) Calculate transitions.
transitions_2P = [
        Rhfs(calc_name = '2p_3',use_ci=True),
        Rbiotransform(use_ci=True,calc_name_initial = '2s_3',calc_name_final = '2p_3', transform_all = True),
        Rtransition(use_ci=True,calc_name_initial = '2s_3',calc_name_final = '2p_3',transition_spec = ['E1'])]

[cmd.execute(workdir = ex1_workdir) for cmd in transitions_2P]
