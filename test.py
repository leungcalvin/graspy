import grasp

ground = grasp.Rcsfgenerate('Xe',
        ['4f(14,*)6s(2,*)','4f(14,*)6s(1,*)6p(1,*)'],
        activeset='8s,8p,7d,6f,6g',
        jlower=0,jhigher=2,exc=0)
ground.execute(workdir='/home/calvin/grasp-wrapper/testdir')
