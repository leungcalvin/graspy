import subprocess
from shutil import copyfile
import os
from os import path

def booltoyesno(boolean):
    if boolean:
        return 'y'
    else:
        return 'n'

def checkSmart(abspath):
    if path.exists(abspath):
        return True
    else:
        root,ext=path.splitext(abspath)
        if path.exists(root+'.out'):
            print(f'Copying {root}.out -> {root}.inp')
            os.rename(root+'.out',root+'.inp')
            return True
        else:
            return False


class Routine(object):
    def __init__(self,name,inputs=[],outputs=[],params=[]):
        self.name = name
        self.inputs=inputs
        self.outputs=outputs
        self.params=params
    def execute(self,workdir):
        self.workdir = workdir
        inputstr = ''.join(p+'\n' for p in self.params)
        for inp in self.inputs:
            assert checkSmart(path.join(workdir,inp)), f'{self.name}: The input file {inp} is missing'
        completedProcess = subprocess.run([self.name],
                input=inputstr,
                #capture_output=False,
                shell=True,
                cwd=workdir,
                check=True,
                encoding='utf8')
        #print(completedProcess.stdout)
        #print(completedProcess.stderr)
        for outp in self.outputs:
            assert path.exists(path.join(workdir,outp)),f'{self.name}: The output file {outp} is missing'


class Rnucleus(Routine):
    def __init__(self,Z,A,neutralMass,I,NDM,NQM):
        params = [str(inp) for inp in [Z,A,'n',neutralMass,I,NDM,NQM]]
        super().__init__(name='rnucleus',
                         inputs=[],
                         outputs=['isodata'],
                         params=params)


coredict = {'He':1,'Ne':2,'Ar':3,'Kr':4,'Xe':5}
class Rcsfgenerate(Routine):
    def __init__(self,core,csflist,activeset,jlower,jhigher,exc):
        # todo: implement non-user specified ordering of orbital iteration
        params = ['u',str(coredict[core])]
        params.extend(csflist)
        params.append('*') # end the CSF list
        params.append(activeset)
        params.append(f'{jlower},{jhigher}')
        params.append(str(exc))
        params.append('n') # don't add another list of CSFs
        super().__init__(name='rcsfgenerate',
                         inputs=['clist.ref'],
                         outputs=['rcsf.out','rcsfgenerate.log'],
                         params=params)
    def execute(self,workdir,writeMR=False):
        # we might need to copy to a multiref file at the end
        super().execute(workdir)
        if writeMR:
            copyfile(path.join(workdir,'rcsf.out'),
                    path.join(workdir,'rcsfmr.inp'))



hamiltoniandict = {'DC':1,'DCB':2}

class Rcsfinteract(Routine):
    def __init__(self,hamiltonian):
        # todo: implement non-user specified ordering of orbital iteration
        params = [str(hamiltoniandict[hamiltonian])]
        super().__init__(name='rcsfinteract',
                         inputs=['rcsfmr.inp','rcsf.inp'],
                         outputs=['rcsf.out'],
                         params=params)

class Rangular(Routine):
    def __init__(self):
        params = ['y']
        super().__init__(name='rangular',
                         inputs=['rcsf.inp'],
                         outputs=['rangular.log'],
                         params=params)
        # Also: rangular will make a bunch of mcp files which we won't keep track of

methoddict = {'Thomas-Fermi':'2','Screened Hydrogenic':'3'}
class Rwfnestimate(Routine):
    def __init__(self,orbdict=None,fallback='Thomas-Fermi'):
        if orbdict == None:
            defaultorbs = {'*':os.path.join(self.workdir,'rwfn.out')}
        params = ['y']
        # make rwfnestimates, but save wildcard entry for last
        for key in orbdict.keys():
            if key != '*':
                params.extend(['1',orbdict[key],key])
        if '*' in orbdict.keys():
            params.extend(['1',orbdict['*'],'*'])

        params.extend([methoddict[fallback],'*']) # after looking up everything possible, make sure all orbitals are estimated according to the fallback method

        super().__init__(name='rwfnestimate',
                         inputs=['isodata','rcsf.inp']+list(orbdict.values()),
                         outputs=['rwfn.inp'],
                         params=params)

weightingmethods = {'Equal':'1','Standard':'5','User [unsupported!]':'9'}
class Rmcdhf(Routine):
    def __init__(self,asfidx,orbs,specorbs,runs,weightingmethod):
        #asf indices should be given as a list whose entries are lists of comma-separated integers which are the serial numbers accepted by GRASP.
        #orbs should be a list of strings like ['6s','6p*'], denoting which orbitals for the RMCDHF routine to optimize
        # runs is an integer denoting the maximum number of SCF iterations.
        params = ['y']
        params.extend( ','.join(str(level) for level in levels) for levels in asfidx)
        params.append(weightingmethods[weightingmethod])
        params.append(','.join(orbs))
        params.append(','.join(specorbs))
        params.append(str(runs))
        super().__init__(name='rmcdhf',
                         inputs = ['isodata','rcsf.inp','rwfn.inp'],
                         outputs = ['rmix.out','rwfn.out','rmcdhf.sum','rmcdhf.log'],
                         params=params)

class Rsave(Routine):
    def __init__(self,calcname):
        super().__init__(name=f'rsave {calcname}',
                         inputs = ['rmix.out','rwfn.out','rmcdhf.sum','rmcdhf.log'],
                         outputs= [f'{calcname}.{ext}' for ext in ['c','m','w','sum','log']],
                         params = [])


class Rci(Routine):
    def __init__(self,calcname,includetransverse,modifyfreq,scalefactor,includevacpol,includenms,includesms,estselfenergy,largestn,asfidx):
        params = ['y',calcname]
        params.extend([booltoyesno(param) for param in [includetransverse,modifyfreq]])
        if includetransverse:
            params.append(str(scalefactor))
        params.extend([booltoyesno(param) for param in [includevacpol,includenms,includesms,estselfenergy]])
        if estselfenergy:
            params.append(str(largestn))
        params.extend( ','.join(str(level) for level in levels) for levels in asfidx)
        super().__init__(name = 'rci',
                         inputs = [f'{calcname}.c',f'{calcname}.w'],
                         outputs=[f'{calcname}.cm',f'{calcname}.csum',f'{calcname}.clog','rci.res'],
                         params = params)


