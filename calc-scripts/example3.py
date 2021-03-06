from graspy.grasp import *

calcdir = '/home/calvinl/graspy/calc-outputs/example3'
initialize(workdir=calcdir)

ref_a = Rcsfgenerate(core='None',ordering = 'Default',
            csflist=['1s(2,i)2s(2,i)2p(3,i)'],
            active_set=[2,2],
            jlower=1,jhigher=5,exc=0,write_csf= 'rcsfmr.inp')
ref_b = Rcsfgenerate(core='None',ordering = 'Default',
            csflist=['1s(2,i)2p(5,i)'],
            active_set = [2,2],
            jlower=1,jhigher=5,exc=0,write_csf= 'rcsfmr.inp')

# The three states which make up the multireference can be added together with arithmetic.
mr = ref_a + ref_b
mr.execute(workdir = calcdir)

MR_DHF =[
        Rnucleus(Z=14,A=28,neutralMass=27.9769271,I=1,NDM=1,NQM=1),
        Rangular(),
        Rwfnestimate(orbdict = None, fallback='Screened Hydrogenic'),
        Rmcdhf([[1,2],[1,2,3,4],[1]],orbs = ['*'],specorbs = ['*'], runs = 100, weighting_method = 'Standard'),
        Rsave('2s22p3_2p5_DF')
        ]
[cmd.execute(workdir = calcdir) for cmd in MR_DHF]
CAS_a = Rcsfgenerate(core='None',ordering = 'Default',
            csflist=['1s(2,*)2s(2,*)2p(3,*)'],
            active_set=[3,3,3],
            jlower=1,jhigher=5,exc=2,write_csf = 'rcsf.inp')

CAS_b = Rcsfgenerate(core='None',ordering = 'Default',
    csflist=['1s(2,*)2p(5,*)'],
    active_set=[3,3,3],
    jlower=1,jhigher=5,exc=2,write_csf = 'rcsf.inp')
CAS_2S_exp = CAS_a + CAS_b
CAS_2S_exp.execute(workdir = calcdir)

interact = Rcsfinteract('Dirac-Coulomb')
interact.execute(workdir=calcdir)

CAS_2S = [
        Rangular(),
        Rwfnestimate(orbdict = {'*':'2s22p3_2p5_DF.w'},fallback = 'Screened Hydrogenic'),
        Rmcdhf([[1,2],[1,2,3,4],[1]],
            orbs = ['3*'],
            specorbs = [' '],
            runs = 100, weighting_method = 'Standard'),
        Rsave('2s22p3_2p5_3')
        ]
[cmd.execute(workdir = calcdir) for cmd in CAS_2S]
indices_try=[[1,2],[1,2,3,4],[1]]

CI_2S = [
        Rci(calc_name='2s22p3_2p5_3',
                include_transverse=True,
                modify_freq=True,
                scale_factor='1.d-6',
                include_vacpol = True,
                include_nms= False,
                include_sms = False,
                est_self_energy= True,
                largest_n = 3,
                asfidx = indices_try),
        JJtoLSJ(calc_name= '2s22p3_2p5_3',use_ci = True, unique = True),
        ]
[cmd.execute(workdir = calcdir) for cmd in CI_2S]

transitions_2P = [
        Rhfs(calc_name = '2s22p3_2p5_3',use_ci=True),
        Rbiotransform(use_ci=True,calc_name_initial = '2s22p3_2p5_3',calc_name_final = '2s22p3_2p5_3', transform_all = True),
        Rtransition(use_ci=True,calc_name_initial = '2s22p3_2p5_3',calc_name_final = '2s22p3_2p5_3',transition_spec = ['E1'])]

[cmd.execute(workdir = calcdir) for cmd in transitions_2P]
