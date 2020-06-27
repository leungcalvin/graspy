from graspy.grasp import *
calc_dir ='/home/calvin/graspy/calc-outputs/yb_basis_set'
initialize(workdir = calc_dir,clist = ['1s',
'2s','2p','3s','3p','3d','4s','4p','4d','5s','5p','4f','5d','6s','6p','5f','6d','7p','8s'])

def gen_nucleus():
    cmds = [Rnucleus(Z=70,A=172,neutralMass=171.936368659,I=0,NDM=0,NQM=0,rms_radius = 5.294,thickness = 2.18)]
    return [cmd.execute(workdir = calc_dir) for cmd in cmds]

def gen_multireference():
    mr = [Rcsfgenerate('Ar','4d(10,c)5s(2,i)5p(6,i)4f(14,i)6s(2,i)',activeset=[6,5,4,4],jlower=0,jhigher=0,exc=0,ordering = 'User specified',write_csf = 'rcsfmr.inp'),
            Rcsfgenerate('Ar','4d(10,c)5s(2,i)5p(6,i)4f(14,i)6s(1,i)6p(1,i)',activeset = [6,5,4,4],jlower=0,jhigher=2,exc=0,ordering = 'User specified',write_csf = 'rcsfmr.inp'),
            Rcsfgenerate('Ar','4d(10,c)5s(2,i)5p(6,i)4f(14,i)6s(1,i)5d(1,i)',activeset = [6,5,4,4],jlower=2,jhigher=4,exc=0,ordering = 'User specified',write_csf = 'rcsfmr.inp'),
            Rcsfgenerate('Ar','4d(10,c)5s(2,i)5p(6,i)4f(14,i)6s(1,i)5f(1,i)',activeset = [6,5,4,5],jlower=4,jhigher=6,exc=0,ordering = 'User specified',write_csf ='rcsfmr.inp'),
]
    combined = mr[0] + mr[1] + mr[2] + mr[3]
    return combined.execute(workdir = calc_dir)
def angular_integration():
    return Rangular().execute(workdir = calc_dir)

def estimate_wavefunctions(grid):
    return Rwfnestimate(grid = grid, orbdict = {'*': '/home/calvin/graspy/calc-scripts/cores/yb_6s2master.w'},fallback = 'Thomas-Fermi').execute(workdir = calc_dir)
master_grid = {'RNT': 2.857142857143E-08, 'H': 0.05, 'HP':0 , 'N': 590}

def run_hartree_fock(grid):
    return Rmcdhf(asfidx = [[1],[1],[1],[1,2],[1,2],[1],[1,2]],orbs = ['*'],specorbs = ['*'], runs = 1000, weighting_method = 'Standard',integration_method = 3,grid = grid).execute(workdir = calc_dir)

gen_nucleus()
gen_multireference()
angular_integration()
estimate_wavefunctions(master_grid)
run_hartree_fock(master_grid)

## construct each level
# stateID="2S12"
# stateIDs+=("$stateID")

#1. Construct MR
# rcsfgenerate << EOF
# u
# 4
# 4d(10,c)5s(2,i)5p(6,i)4f(14,i)6s(1,i)
# 4d(10,c)5s(2,i)5p(6,i)4f(14,i)6p(1,i)
# 4d(10,c)5s(2,i)5p(6,i)4f(14,i)5d(1,i)
# 4d(10,c)5s(2,i)5p(6,i)4f(13,i)6s(2,i)
# 4d(10,c)5s(2,i)5p(6,i)4f(13,i)6s(1,i)5d(1,i)
# 4d(10,c)5s(2,i)5p(6,i)4f(13,i)6p(2,i)
# 4d(10,c)5s(2,i)5p(6,i)4f(13,i)5d(2,i)
#
# 6s,6p,5d,4f
# 1,1
# 0
# n
# EOF
# # cp rcsf.out rcsfmr.inp
# cp rcsf.out rcsf.inp

# Leave only target levels
# rsave $stateID
# printf "$stateID\ny" | rasfsplit
# list=("even1")
# printf "" > rcsf.inp # initialize file
# sed -n '1,5p' < "$stateID.c" >> rcsf.inp
# sed -n "6,$(sed -n '$=' ""${stateID}_${list[0]}.c"")p" < "${stateID}_${list[0]}.c" >> rcsf.inp
# for f in ${list[@]:1}; do
# 	printf " *\n" >> rcsf.inp
# 	sed -n "6,$(sed -n '$=' ""${stateID}_$f.c"")p" < "${stateID}_$f.c" >> rcsf.inp
# done
# rm "${stateID}"*.c
# # cp rcsf.inp rcsfmr.inp
# cp rcsf.inp rcsfmr.inp.$stateID
#
# #2. Excite MR
# rcsfgenerate << EOF
# u
# 4
# 4d(10,c)5s(2,i)5p(6,i)4f(14,13)6s(1,*)
# 4d(10,c)5s(2,i)5p(6,i)4f(14,13)6p(1,*)
# 4d(10,c)5s(2,i)5p(6,i)4f(14,13)5d(1,*)
# 4d(10,c)5s(2,i)5p(6,i)4f(13,12)6s(2,*)
# 4d(10,c)5s(2,i)5p(6,i)4f(13,12)6s(1,*)5d(1,*)
# 4d(10,c)5s(2,i)5p(6,i)4f(13,12)6p(2,*)
# 4d(10,c)5s(2,i)5p(6,i)4f(13,12)5d(2,*)
#
# 7s,7p,7d,7f
# 1,1
# 2
# n
# EOF
# cp rcsf.out rcsf.inp
#
# # Leave only target levels
# rsave $stateID
# printf "$stateID\ny" | rasfsplit
# # list=("even1")
# printf "" > rcsf.inp # initialize file
# sed -n '1,5p' < "$stateID.c" >> rcsf.inp
# sed -n "6,$(sed -n '$=' ""${stateID}_${list[0]}.c"")p" < "${stateID}_${list[0]}.c" >> rcsf.inp
# for f in ${list[@]:1}; do
# 	printf " *\n" >> rcsf.inp
# 	sed -n "6,$(sed -n '$=' ""${stateID}_$f.c"")p" < "${stateID}_$f.c" >> rcsf.inp
# done
# rm "${stateID}"*.c
# # cp rcsf.inp rcsf.inp.$stateID
# cp rcsf.inp rcsf.inp.$stateID
#
#
# stateID="2Dduplet"
# stateIDs+=("$stateID")
#
# #1. Construct MR
# rcsfgenerate << EOF
# u
# 4
# 4d(10,c)5s(2,i)5p(6,i)4f(14,i)6s(1,i)
# 4d(10,c)5s(2,i)5p(6,i)4f(14,i)6p(1,i)
# 4d(10,c)5s(2,i)5p(6,i)4f(14,i)5d(1,i)
# 4d(10,c)5s(2,i)5p(6,i)4f(13,i)6s(2,i)
# 4d(10,c)5s(2,i)5p(6,i)4f(13,i)6s(1,i)5d(1,i)
# 4d(10,c)5s(2,i)5p(6,i)4f(13,i)6p(2,i)
# 4d(10,c)5s(2,i)5p(6,i)4f(13,i)5d(2,i)
#
# 6s,6p,5d,4f
# 3,5
# 0
# n
# EOF
# # cp rcsf.out rcsfmr.inp
# cp rcsf.out rcsf.inp
#
# # Leave only target levels
# rsave $stateID
# printf "$stateID\ny" | rasfsplit
# list=("even2" "even1")
# printf "" > rcsf.inp # initialize file
# sed -n '1,5p' < "$stateID.c" >> rcsf.inp
# sed -n "6,$(sed -n '$=' ""${stateID}_${list[0]}.c"")p" < "${stateID}_${list[0]}.c" >> rcsf.inp
# for f in ${list[@]:1}; do
# 	printf " *\n" >> rcsf.inp
# 	sed -n "6,$(sed -n '$=' ""${stateID}_$f.c"")p" < "${stateID}_$f.c" >> rcsf.inp
# done
# rm "${stateID}"*.c
# # cp rcsf.inp rcsfmr.inp
# cp rcsf.inp rcsfmr.inp.$stateID
#
# #2. Excite MR
# rcsfgenerate << EOF
# u
# 4
# 4d(10,c)5s(2,i)5p(6,i)4f(14,13)6s(1,*)
# 4d(10,c)5s(2,i)5p(6,i)4f(14,13)6p(1,*)
# 4d(10,c)5s(2,i)5p(6,i)4f(14,13)5d(1,*)
# 4d(10,c)5s(2,i)5p(6,i)4f(13,12)6s(2,*)
# 4d(10,c)5s(2,i)5p(6,i)4f(13,12)6s(1,*)5d(1,*)
# 4d(10,c)5s(2,i)5p(6,i)4f(13,12)6p(2,*)
# 4d(10,c)5s(2,i)5p(6,i)4f(13,12)5d(2,*)
#
# 7s,7p,7d,7f
# 3,5
# 2
# n
# EOF
# cp rcsf.out rcsf.inp
#
# # Leave only target levels
# rsave $stateID
# printf "$stateID\ny" | rasfsplit
# # list=("even1")
# printf "" > rcsf.inp # initialize file
# sed -n '1,5p' < "$stateID.c" >> rcsf.inp
# sed -n "6,$(sed -n '$=' ""${stateID}_${list[0]}.c"")p" < "${stateID}_${list[0]}.c" >> rcsf.inp
# for f in ${list[@]:1}; do
# 	printf " *\n" >> rcsf.inp
# 	sed -n "6,$(sed -n '$=' ""${stateID}_$f.c"")p" < "${stateID}_$f.c" >> rcsf.inp
# done
# rm "${stateID}"*.c
# cp rcsf.inp rcsf.inp.$stateID
#
#
# stateID="2F72"
# stateIDs+=("$stateID")
#
# #1. Construct MR
# rcsfgenerate << EOF
# u
# 4
# 4d(10,c)5s(2,i)5p(6,i)4f(14,i)6s(1,i)
# 4d(10,c)5s(2,i)5p(6,i)4f(14,i)6p(1,i)
# 4d(10,c)5s(2,i)5p(6,i)4f(14,i)5d(1,i)
# 4d(10,c)5s(2,i)5p(6,i)4f(13,i)6s(2,i)
# 4d(10,c)5s(2,i)5p(6,i)4f(13,i)6s(1,i)5d(1,i)
# 4d(10,c)5s(2,i)5p(6,i)4f(13,i)6p(2,i)
# 4d(10,c)5s(2,i)5p(6,i)4f(13,i)5d(2,i)
#
# 6s,6p,5d,4f
# 7,7
# 0
# n
# EOF
# # cp rcsf.out rcsfmr.inp
# cp rcsf.out rcsfmr.inp.$stateID
# # cp rcsf.out rcsf.inp
#
# # # Leave only target levels
# # rsave $stateID
# # printf "$stateID\ny" | rasfsplit
# # list=("even1")
# # printf "" > rcsf.inp # initialize file
# # sed -n '1,5p' < "$stateID.c" >> rcsf.inp
# # sed -n "6,$(sed -n '$=' ""${stateID}_${list[0]}.c"")p" < "${stateID}_${list[0]}.c" >> rcsf.inp
# # for f in ${list[@]:1}; do
# # 	printf " *\n" >> rcsf.inp
# # 	sed -n "6,$(sed -n '$=' ""${stateID}_$f.c"")p" < "${stateID}_$f.c" >> rcsf.inp
# # done
# # rm "${stateID}"*.c
# # cp rcsf.inp rcsfmr.inp
#
# #2. Excite MR
# rcsfgenerate << EOF
# u
# 4
# 4d(10,c)5s(2,i)5p(6,i)4f(14,13)6s(1,*)
# 4d(10,c)5s(2,i)5p(6,i)4f(14,13)6p(1,*)
# 4d(10,c)5s(2,i)5p(6,i)4f(14,13)5d(1,*)
# 4d(10,c)5s(2,i)5p(6,i)4f(13,12)6s(2,*)
# 4d(10,c)5s(2,i)5p(6,i)4f(13,12)6s(1,*)5d(1,*)
# 4d(10,c)5s(2,i)5p(6,i)4f(13,12)6p(2,*)
# 4d(10,c)5s(2,i)5p(6,i)4f(13,12)5d(2,*)
#
# 7s,7p,7d,7f
# 7,7
# 2
# n
# EOF
# cp rcsf.out rcsf.inp
#
# # Leave only target levels
# rsave $stateID
# printf "$stateID\ny" | rasfsplit
# list=("odd1")
# printf "" > rcsf.inp # initialize file
# sed -n '1,5p' < "$stateID.c" >> rcsf.inp
# sed -n "6,$(sed -n '$=' ""${stateID}_${list[0]}.c"")p" < "${stateID}_${list[0]}.c" >> rcsf.inp
# for f in ${list[@]:1}; do
# 	printf " *\n" >> rcsf.inp
# 	sed -n "6,$(sed -n '$=' ""${stateID}_$f.c"")p" < "${stateID}_$f.c" >> rcsf.inp
# done
# rm "${stateID}"*.c
# cp rcsf.inp rcsf.inp.$stateID
#
#
#
# ## Combine rcsf from each level
# cat rcsfmr.inp.${stateIDs[0]} > rcsfmr.inp
# for f in ${stateIDs[@]:1}; do
# 	printf " *\n" >> rcsfmr.inp
# 	sed -n "6,$(sed -n '$=' ""rcsfmr.inp.$f"")p" < "rcsfmr.inp.$f" >> rcsfmr.inp
# done
#
# cat rcsf.inp.${stateIDs[0]} > rcsf.inp
# for f in ${stateIDs[@]:1}; do
# 	printf " *\n" >> rcsf.inp
# 	sed -n "6,$(sed -n '$=' ""rcsf.inp.$f"")p" < "rcsf.inp.$f" >> rcsf.inp
# done
#
#
#
# ## run the rest part
# stateID=$stateID_all
#
# printf "1" | rcsfinteract
# cp rcsf.out rcsf.inp
# # printf 'y' | rangular
# # MPI_NP=4 # # of processors for MPI
# printf 'y' | mpirun -np $MPI_NP rangular_mpi
# # rwfnestimate_input="input_${stateID}_rwfnestimate"
# # rwfnestimate < $rwfnestimate_input
# # cat $rwfnestimate_input
# cp "rwfn_8spdf.out" rwfn.inp
# # rmcdhf < $rmcdhf_input
# sleep 2s
# mpirun -np $MPI_NP rmcdhf_mpi << EOF
# y
# 1
# 1
# 1
# 1
# 5
#
# *
# 1000
# EOF
#
# rsave $stateID
#
# printf "$stateID\nn\ny\ny" | jj2lsj
# printf "$stateID\nn\n0.01\ny" | rmixextract >& $stateID.rmix # DHF
#
# rlevelseV $stateID.m |& tee $stateID.rlevels.DHF
#
# sleep 2s
# ## CI
# # rci << EOF
# # MPI_NP=4 # # of processors for MPI
# mpirun -np $MPI_NP rci_mpi << EOF
# y
# $stateID
# y
# y
# 1.d-6
# y
# n
# n
# y
# 8
# 1
# 1
# 1
# 1
# EOF
#
# printf "$stateID\ny\ny\ny" | jj2lsj
# printf "$stateID\ny\n0.01\ny" | rmixextract >& $stateID.crmix # CI
#
# rlevelseV $stateID.cm |& tee $stateID.rlevels.CI

# printf "y\n$stateID\ny\ny\nn" | ris4
# printf "y\n$stateID\ny" | redf


