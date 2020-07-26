from graspy.grasp import *

def gen_nucleus(calc_dir,Z):
    cmds = [Rnucleus(Z=Z,A=172,neutralMass=171.936368659,I=0,NDM=0,NQM=0,rms_radius = 5.294,thickness = 2.18)]
    return [cmd.execute(workdir = calc_dir) for cmd in cmds]

def angular_integration(calc_dir):
    return Rangular().execute(workdir = calc_dir)

def estimate_wavefunctions(calc_dir,previous_rwfn, grid):
    return Rwfnestimate(grid = grid, orbdict = {'*': previous_rwfn},fallback = 'Thomas-Fermi').execute(workdir = calc_dir)

master_grid = {'RNT': 2.857142857143E-08, 'H': 0.05, 'HP':0 , 'N': 590}

def run_hartree_fock_save(calc_dir,grid,orbs = ['*'],specorbs = ['*'],integration_method = None,calc_name = None):
    #out = Rmcdhf(asfidx = [[1],[1],[1],[1,2],[1,2],[1,2],[1],[1,2],[1]],orbs = orbs, specorbs = specorbs, runs = 20000, weighting_method = 'Standard',integration_method = integration_method,grid = grid).execute(workdir = calc_dir)
    out = Rmcdhf(asfidx = [[1],[1],[1],[1,2],[1,2],[1,2,3,4,5],[1],[1,2,3,4],[1,2,3,4],[1,2,3]],orbs = orbs, specorbs = specorbs, runs = 20000, weighting_method = 'Standard',integration_method = integration_method,grid = grid).execute(workdir = calc_dir)

    if calc_name is not None:
        Rsave(calc_name).execute(workdir = calc_dir)
        return os.path.join(calc_dir,calc_name) + '.w'
    else:
        return out
def run_ci(calc_dir,calc_name,largest_n):
    return Rci(calcname = calc_name,
               includetransverse=True,
               modifyfreq=True,
               scalefactor='1.d-6',
               includevacpol = True,
               includenms= False,
               includesms = False,
               estselfenergy= True,
               largestn = largest_n,
               asfidx = [[1],[1],[1],[1,2],[1,2],[1,2,3,4,5],[1],[1,2,3,4],[1,2,3,4],[1,2,3]]).execute(workdir = calc_dir)

def csfs_open_maxn(calc_dir,write_csf,maxn,exc):
    # max_n: a list of 4 quantum numbers representing the maximum s,p,d, and f shells.
    # 4f14 6s2  with 6s excitable
    # 4f14 6s 6p with 6p excitable
    # 4f14 6s 5d with 5d excitable
    # 4f13 6s2 5d with 4f excitable
    mr =   [Rcsfgenerate('Kr','4d(10,c)5s(2,i)5p(6,i)4f(14,i)6s(2,*)',
            activeset=[maxn[0],5,4,4],jlower=0,jhigher=2,exc=exc,
            ordering = 'User specified',write_csf = write_csf),

            Rcsfgenerate('Kr','4d(10,c)5s(2,i)5p(6,i)4f(14,i)6s(1,i)6p(1,*)',         activeset = [6,maxn[1],5,4],jlower=0,jhigher=4,exc=exc,
            ordering = 'User specified',write_csf = write_csf),

            Rcsfgenerate('Kr','4d(10,c)5s(2,i)5p(6,i)4f(14,i)6s(1,i)5d(1,*)',         activeset = [6,5,maxn[2],4],jlower=2,jhigher=6,exc=exc,
            ordering = 'User specified',write_csf = write_csf),

            Rcsfgenerate('Kr','4d(10,c)5s(2,i)5p(6,i)4f(13,*)6s(2,i)5d(1,i)',         activeset = [6,5,5,maxn[3]],jlower=4,jhigher=10,exc=exc,
            ordering = 'User specified',write_csf = write_csf)]
    combined = mr[0] + mr[1] + mr[2] + mr[3]
    return combined.execute(workdir = calc_dir)

def csfs_open_maxn_safronova(calc_dir,write_csf,maxn,exc):
    # max_n: a list of 4 quantum numbers representing the maximum s,p,d, and f shells.
    # 4f14 6s2  with 6s excitable
    # 4f14 6s 6p with 6p excitable
    # 4f14 6s 5d with 5d excitable
    # 4f13 6s2 5d with 4f excitable
    # 4f13 6s2 6p with 6p,4f excitable
    mr =   [Rcsfgenerate('Kr','4d(10,c)5s(2,i)5p(6,i)4f(14,i)6s(2,*)',
            activeset=[maxn[0],5,4,4],jlower=0,jhigher=2,exc=exc,
            ordering = 'User specified',write_csf = write_csf),

            Rcsfgenerate('Kr','4d(10,c)5s(2,i)5p(6,i)4f(14,i)6s(1,i)6p(1,*)',         activeset = [6,maxn[1],5,4],jlower=0,jhigher=4,exc=exc,
            ordering = 'User specified',write_csf = write_csf),

            Rcsfgenerate('Kr','4d(10,c)5s(2,i)5p(6,i)4f(14,i)6s(1,i)5d(1,*)',         activeset = [6,5,maxn[2],4],jlower=2,jhigher=6,exc=exc,
            ordering = 'User specified',write_csf = write_csf),

            Rcsfgenerate('Kr','4d(10,c)5s(2,i)5p(6,i)4f(13,*)6s(2,i)5d(1,i)',         activeset = [6,5,5,maxn[3]],jlower=4,jhigher=10,exc=exc,
            ordering = 'User specified',write_csf = write_csf),
            Rcsfgenerate('Kr','4d(10,c)5s(2,i)5p(6,i)4f(13,*)6s(2,i)6p(1,*)',         activeset = [6,maxn[1],4,maxn[3]],jlower=0,jhigher=4,exc=exc,
            ordering = 'User specified',write_csf = write_csf)]
    combined = mr[0] + mr[1] + mr[2] + mr[3] + mr[4]
    return combined.execute(workdir = calc_dir)

calc_dirs =['/home/calvin/graspy/calc-outputs/yb_basis_set/spec',
            '/home/calvin/graspy/calc-outputs/yb_basis_set/6spd5f',
            '/home/calvin/graspy/calc-outputs/yb_basis_set/7spd6f',
            '/home/calvin/graspy/calc-outputs/yb_basis_set/8spdf',]

rwfn_out_files = ['/home/calvin/graspy/calc-scripts/cores/yb_6s2master.w']
################
# Spectroscopic Orbitals
################
initialize(workdir = calc_dirs[0],clist = ['1s','2s','2p','3s','3p','3d','4s','4p','4d','5s','5p','4f','5d','6s','6p','5f','6d','7p','8s'])
gen_nucleus(calc_dirs[0],Z = 70)
csfs_open_maxn_safronova(calc_dirs[0],maxn = [6,6,5,4],exc = 0,write_csf = 'rcsfmr.inp')
angular_integration(calc_dirs[0])
estimate_wavefunctions(calc_dirs[0],previous_rwfn = rwfn_out_files[-1], grid = master_grid) # use yb_6s2master.w
rwfn_out_files.append(run_hartree_fock_save(calc_dirs[0],master_grid,orbs = ['1*','2*','3*','4s*','4p*','4d*'],specorbs = ['*'],integration_method = None,calc_name = 'spec'))
# input('1-4d OK? ')
estimate_wavefunctions(calc_dirs[0],previous_rwfn = rwfn_out_files[-1], grid = master_grid)
rwfn_out_files.append(run_hartree_fock_save(calc_dirs[0],master_grid,orbs = ['1*','2*','3*','4s*','4p*','4d*','5s*','5p*'],specorbs = ['*'],integration_method = 4,calc_name = 'spec')) # method 4 works.
# input('5sp OK? ')
estimate_wavefunctions(calc_dirs[0],previous_rwfn = rwfn_out_files[-1], grid = master_grid)
rwfn_out_files.append(run_hartree_fock_save(calc_dirs[0],master_grid,orbs = ['1*','2*','3*','4s*','4p*','4d*','5s*','5p*','6s*'],specorbs = ['*'],integration_method = 4,calc_name = 'spec'))
# input('6s OK?')
estimate_wavefunctions(calc_dirs[0],previous_rwfn = rwfn_out_files[-1], grid = master_grid)
rwfn_out_files.append(run_hartree_fock_save(calc_dirs[0],master_grid,orbs = ['1*','2*','3*','4s*','4p*','4d*','5s*','5p*','6s*','4f*'],specorbs = ['*'],integration_method = 4,calc_name = 'spec'))
# input('4f OK?')
estimate_wavefunctions(calc_dirs[0],previous_rwfn = rwfn_out_files[-1], grid = master_grid)
rwfn_out_files.append(run_hartree_fock_save(calc_dirs[0],master_grid,orbs = ['1*','2*','3*','4s*','4p*','4d*','5s*','5p*','6s*','4f*','5d*'],specorbs = ['*'],integration_method = 4,calc_name = 'spec'))
# input('5d OK?')
estimate_wavefunctions(calc_dirs[0],previous_rwfn = rwfn_out_files[-1], grid = master_grid)
rwfn_out_files.append(run_hartree_fock_save(calc_dirs[0],master_grid,orbs = ['*'],specorbs = ['*'],integration_method = 4,calc_name = 'spec'))
# input('All spectroscopic orbitals OK?') they all look fine!
#################
# Add 6d and 5f
#################
initialize(workdir = calc_dirs[1],clist = ['1s','2s','2p','3s','3p','3d','4s','4p','4d','5s','5p','4f','5d','6s','6p','5f','6d','7p','8s'])
csfs_open_maxn_safronova(calc_dirs[1],maxn = [6,6,6,5],exc = 1,write_csf = 'rcsf.inp')
angular_integration(calc_dirs[1])
z_iter_rwfn = []
import numpy as np
for Z in np.linspace(71,70,num=5):
    gen_nucleus(calc_dirs[1],Z = Z)
    estimate_wavefunctions(calc_dirs[1],previous_rwfn = rwfn_out_files[-1], grid = master_grid)
    z_iter_rwfn.append(run_hartree_fock_save(calc_dirs[1],master_grid,orbs = ['6d*'],specorbs =[] ,integration_method = 3,calc_name = '6spd5f'))
    estimate_wavefunctions(calc_dirs[1],previous_rwfn = z_iter_rwfn[-1], grid = master_grid)
    z_iter_rwfn.append(run_hartree_fock_save(calc_dirs[1],master_grid,orbs = ['5f*'],specorbs = [],integration_method = 3,calc_name = '6spd5f'))
estimate_wavefunctions(calc_dirs[1],previous_rwfn = rwfn_out_files[-1], grid = master_grid)
# input(f'6d 5f OK Z = {Z}?') both look fine, a bit close to continuum state.
rwfn_out_files.append(z_iter_rwfn[-1])
run_ci(calc_dir = calc_dirs[1],calc_name = '6spd5f',largest_n = 7)
###############
# Add 7s
# Joon suggests: 7s2 and 6s2 in configuration.
###############
initialize(workdir = calc_dirs[2],clist = ['1s','2s','2p','3s','3p','3d','4s','4p','4d','5s','5p','4f','5d','6s','6p','5f','6d','7p','8s'])
#csfs_open_maxn_safronova(calc_dirs[2],maxn = [7,7,7,6],exc = 1,write_csf = 'rcsf.inp')
Rcsfgenerate('Kr','4d(10,c)5s(2,i)5p(6,i)4f(14,i)6s(2,*)',
        activeset=[7,5,4,4],jlower=0,jhigher=2,exc=1,
        ordering = 'User specified',write_csf = 'rcsf.inp').execute(workdir = calc_dirs[2])
angular_integration(calc_dirs[2])
specorbs = ['1*','2*','3*','4s*','4p*','4d*','5s*','5p*','6s*','4f*','5d*']
gen_nucleus(calc_dirs[2],Z = Z)
estimate_wavefunctions(calc_dirs[2],previous_rwfn = rwfn_out_files[-1], grid = master_grid) # use result from 6s6p5d5f
rwfn_out_files.append(run_hartree_fock_save(calc_dirs[2],master_grid,orbs = ['7s*'],specorbs = [],integration_method = None,calc_name = '7s6pd5f'))


