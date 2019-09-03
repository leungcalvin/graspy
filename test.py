import grasp

nucleus= grasp.Rnucleus(Z=70,A=172,neutralMass=171.936378,I=0,NDM=0,NQM=0)
ground = grasp.Rcsfgenerate('Xe',
        ['4f(14,*)6s(2,*)','4f(14,*)6s(1,*)6p(1,*)'],
        activeset='8s,8p,7d,6f,6g',
        jlower=0,jhigher=4,exc=0)
angular= grasp.Rangular()
ansatz = grasp.Rwfnestimate(orbdict={'*':'/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/6s2/master.w'},fallback='Thomas-Fermi')
asfindices = [[1],[1],[1,2],[1]]
scf = grasp.Rmcdhf(asfindices,orbs=['6p'],specorbs=['6p'],weightingmethod='Standard',runs=1000)
savecalc = grasp.Rsave('both')
rcicalc = grasp.Rci(calcname='both',
        includetransverse=True,
        modifyfreq=True,
        scalefactor='1.d-6',
        includevacpol = True,
        includenms= True,
        includesms = True,
        estselfenergy= True,
        largestn = 8,
        asfidx = asfindices)

testdir = '/home/calvin/grasp-wrapper/testdir'
nucleus.execute(workdir=testdir)
ground.execute( workdir=testdir,writeMR=True)
angular.execute(workdir=testdir)
ansatz.execute( workdir=testdir)
scf.execute( workdir =testdir)
savecalc.execute( workdir = testdir)


rcicalc.execute(workdir = testdir)
