from graspy.grasp import *

def gen_nucleus(calc_dir,Z):
    cmds = [Rnucleus(Z=Z,A=172,neutralMass=171.936368659,I=0,NDM=0,NQM=0,rms_radius = 5.294,thickness = 2.18)]
    return [cmd.execute(workdir = calc_dir) for cmd in cmds]

def angular_integration(calc_dir):
    return Rangular().execute(workdir = calc_dir)

def csfs_6s6p5d5f_mr(calc_dir,write_csf):
    # 4f14 6s1 + [6s,6p,5d,5f], no excitations.
    mr =   [Rcsfgenerate('Ar','4d(10,c)5s(2,i)5p(6,i)4f(14,i)6s(2,i)',activeset=[6,5,4,4],jlower=0,jhigher=0,exc=0,ordering = 'User specified',write_csf = write_csf),
            Rcsfgenerate('Ar','4d(10,c)5s(2,i)5p(6,i)4f(14,i)6s(1,i)6p(1,i)',activeset = [6,6,4,4],jlower=0,jhigher=4,exc=0,ordering = 'User specified',write_csf = write_csf),
            Rcsfgenerate('Ar','4d(10,c)5s(2,i)5p(6,i)4f(14,i)6s(1,i)5d(1,i)',activeset = [6,5,5,4],jlower=2,jhigher=6,exc=0,ordering = 'User specified',write_csf = write_csf),
            Rcsfgenerate('Ar','4d(10,c)5s(2,i)5p(6,i)4f(14,i)6s(1,i)5f(1,i)',activeset = [6,5,4,5],jlower=4,jhigher=8,exc=0,ordering = 'User specified',write_csf =write_csf),
    ]
    combined = mr[0] + mr[1] + mr[2] + mr[3]
    return combined.execute(workdir = calc_dir)

def csfs_6spdf(calc_dir,write_csf):
    # 4f14 6s1 + [6s,6p,6d,6f], no excitations.
    mr =   [Rcsfgenerate('Ar','4d(10,c)5s(2,i)5p(6,i)4f(14,i)6s(2,i)',activeset=[6,5,4,4],jlower=0,jhigher=0,exc=0,ordering = 'User specified',write_csf = write_csf),
            Rcsfgenerate('Ar','4d(10,c)5s(2,i)5p(6,i)4f(14,i)6s(1,i)6p(1,i)',activeset = [6,6,4,4],jlower=0,jhigher=4,exc=0,ordering = 'User specified',write_csf = write_csf),
            Rcsfgenerate('Ar','4d(10,c)5s(2,i)5p(6,i)4f(14,i)6s(1,i)5d(1,*)',activeset = [6,5,6,4],jlower=2,jhigher=6,exc=1,ordering = 'User specified',write_csf = write_csf),
            Rcsfgenerate('Ar','4d(10,c)5s(2,i)5p(6,i)4f(14,i)6s(1,i)5f(1,*)',activeset = [6,5,4,6],jlower=4,jhigher=8,exc=1,ordering = 'User specified',write_csf =write_csf),
    ]
    combined = mr[0] + mr[1] + mr[2] + mr[3]
    return combined.execute(workdir = calc_dir)

def csfs_7spdf(calc_dir,write_csf):
    # 4f14 6s1 + [6s,6p,5d,5f] base configurations. for each config: allow promotion of one electron to n+2.
    mr =   [Rcsfgenerate('Ar','4d(10,c)5s(2,i)5p(6,i)4f(14,i)6s(2,*)',
            activeset=[7,5,4,4],jlower=0,jhigher=2,exc=1,
            ordering = 'User specified',write_csf = write_csf),

            Rcsfgenerate('Ar','4d(10,c)5s(2,i)5p(6,i)4f(14,i)6s(1,i)6p(1,*)',         activeset = [6,7,4,4],jlower=0,jhigher=4,exc=1,
            ordering = 'User specified',write_csf = write_csf),

            Rcsfgenerate('Ar','4d(10,c)5s(2,i)5p(6,i)4f(14,i)6s(1,i)5d(1,*)',         activeset = [6,5,7,4],jlower=2,jhigher=6,exc=1,
            ordering = 'User specified',write_csf = write_csf),

            Rcsfgenerate('Ar','4d(10,c)5s(2,i)5p(6,i)4f(14,i)6s(1,i)5f(1,*)',         activeset = [6,5,4,7],jlower=4,jhigher=8,exc=1,
            ordering = 'User specified',write_csf = write_csf)]
    combined = mr[0] + mr[1] + mr[2] + mr[3]
    return combined.execute(workdir = calc_dir)

def csfs_8spdf(calc_dir,write_csf):
    # 4f14 6s1 + [6s,6p,5d,5f] base configurations. for each config: allow promotion of one electron to n+2.
    mr =   [Rcsfgenerate('Ar','4d(10,c)5s(2,i)5p(6,i)4f(14,i)6s(2,*)',
            activeset=[8,5,4,4],jlower=0,jhigher=2,exc=1,
            ordering = 'User specified',write_csf = write_csf),

            Rcsfgenerate('Ar','4d(10,c)5s(2,i)5p(6,i)4f(14,i)6s(1,i)6p(1,*)',         activeset = [6,8,4,4],jlower=0,jhigher=4,exc=1,
            ordering = 'User specified',write_csf = write_csf),

            Rcsfgenerate('Ar','4d(10,c)5s(2,i)5p(6,i)4f(14,i)6s(1,i)5d(1,*)',         activeset = [6,5,8,4],jlower=2,jhigher=6,exc=1,
            ordering = 'User specified',write_csf = write_csf),

            Rcsfgenerate('Ar','4d(10,c)5s(2,i)5p(6,i)4f(14,i)6s(1,i)5f(1,*)',         activeset = [6,5,4,8],jlower=4,jhigher=8,exc=1,
            ordering = 'User specified',write_csf = write_csf)]
    combined = mr[0] + mr[1] + mr[2] + mr[3]
    return combined.execute(workdir = calc_dir)

def estimate_wavefunctions(calc_dir,previous_rwfn, grid):
    return Rwfnestimate(grid = grid, orbdict = {'*': previous_rwfn},fallback = 'Thomas-Fermi').execute(workdir = calc_dir)

master_grid = {'RNT': 2.857142857143E-08, 'H': 0.05, 'HP':0 , 'N': 590}

def run_hartree_fock_save(calc_dir,grid,orbs = ['*'],specorbs = ['*'],integration_method = None,calc_name = None):
    out = Rmcdhf(asfidx = [[1],[1],[1],[1,2],[1,2],[1,2],[1],[1,2],[1]],orbs = orbs, specorbs = specorbs, runs = 1000, weighting_method = 'Standard',integration_method = integration_method,grid = grid).execute(workdir = calc_dir)

    if calc_name is not None:
        Rsave(calc_name).execute(workdir = calc_dir)
        return os.path.join(calc_dir,calc_name) + '.w'
    else:
        return out

calc_dirs =['/home/calvin/graspy/calc-outputs/yb_basis_set/6s6p5d5f',
            '/home/calvin/graspy/calc-outputs/yb_basis_set/6spdf',
            '/home/calvin/graspy/calc-outputs/yb_basis_set/7spdf',
            '/home/calvin/graspy/calc-outputs/yb_basis_set/8spdf',]

rwfn_out_files = ['/home/calvin/graspy/calc-scripts/cores/yb_6s2master.w']
initialize(workdir = calc_dirs[0],clist = ['1s','2s','2p','3s','3p','3d','4s','4p','4d','5s','5p','4f','5d','6s','6p','5f','6d','7p','8s'])
gen_nucleus(calc_dirs[0],Z = 70)
csfs_6s6p5d5f_mr(calc_dirs[0],write_csf = 'rcsfmr.inp')
angular_integration(calc_dirs[0])
estimate_wavefunctions(calc_dirs[0],previous_rwfn = rwfn_out_files[-1], grid = master_grid) # use yb_6s2master.w
core_file = run_hartree_fock_save(calc_dirs[0],master_grid,integration_method = 4,calc_name = 'core') # methods 3 and 4 work!
rwfn_out_files.append(core_file)

initialize(workdir = calc_dirs[1],clist = ['1s','2s','2p','3s','3p','3d','4s','4p','4d','5s','5p','4f','5d','6s','6p','5f','6d','7p','8s'])
gen_nucleus(calc_dirs[1],Z = 70)
csfs_6spdf(calc_dirs[1],write_csf = 'rcsf.inp')
angular_integration(calc_dirs[1])
estimate_wavefunctions(calc_dirs[1],previous_rwfn = rwfn_out_files[-1], grid = master_grid) # use result from 6s6p5d5f
rwfn_out_files.append(run_hartree_fock_save(calc_dirs[1],master_grid,orbs = ['6d*','6f*'],specorbs = [''],integration_method = None,calc_name = '6spdf'))

initialize(workdir = calc_dirs[2],clist = ['1s','2s','2p','3s','3p','3d','4s','4p','4d','5s','5p','4f','5d','6s','6p','5f','6d','7p','8s'])
gen_nucleus(calc_dirs[2],Z = 70)
csfs_7spdf(calc_dirs[2],write_csf = 'rcsf.inp')
angular_integration(calc_dirs[2])
estimate_wavefunctions(calc_dirs[2],previous_rwfn = rwfn_out_files[-1], grid = master_grid) # use result from 6spdf
rwfn_out_files.append(run_hartree_fock_save(calc_dirs[2],master_grid,orbs = ['7s*','7p*','7d*','7f*'],specorbs = [''],integration_method = None,calc_name = '7spdf'))

initialize(workdir = calc_dirs[3],clist = ['1s','2s','2p','3s','3p','3d','4s','4p','4d','5s','5p','4f','5d','6s','6p','5f','6d','7p','8s'])
gen_nucleus(calc_dirs[3],Z = 70)
csfs_8spdf(calc_dirs[3],write_csf = 'rcsf.inp')
angular_integration(calc_dirs[3])
estimate_wavefunctions(calc_dirs[3],previous_rwfn = rwfn_out_files[-1], grid = master_grid) # use result from 6spdf
rwfn_out_files.append(run_hartree_fock_save(calc_dirs[3],master_grid,orbs = ['8s*','8p*','8d*','8f*'],specorbs = [''],integration_method = None,calc_name = '8spdf'))
print(rwfn_out_files)
