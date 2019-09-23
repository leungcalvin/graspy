import pdb
import numpy as np
import pandas as pd
import subprocess
from shutil import copyfile,move
import os
from os import path

def initialize(workdir,clist):
    """
    Makes working directory and populates it with an orbital ordering by writing a workdir/clist.ref file
    Inputs:
    -------
    workdir (str): Absolute path to new directory
    clist (list of str): e.g. ['1s','2s','2p']
    -------
    """
    if not os.path.exists(workdir):
        os.makedirs(workdir)
    with open(os.path.join(workdir,'clist.ref'),'w+') as clistfile:
        clistfile.write(str(''.join(string+'\n') for string in clist))

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
        print(f'{name}: {params}')
        print(f'{name}: {inputs}')
        print(f'{name}: {outputs}')
        self.params=params
        self.hash=hash(self)

    def execute(self,workdir):
        self.workdir = workdir
        inputstr = ''.join(p+'\n' for p in self.params)
        for inp in self.inputs:
            assert checkSmart(path.join(workdir,inp)), f'{self.name}: The input file {inp} is missing'
        process = subprocess.run([f'{self.name} | tee readout.temp'], # capture output to a temp file
                input=inputstr,
                capture_output=False,
                #stdout=subprocess.PIPE,
                shell=True,
                cwd=workdir,
                check=True,
                encoding='utf8')
        for outp in self.outputs:
            assert path.exists(path.join(workdir,outp)),f'{self.name}: The output file {outp} is missing'
        # remove temp file
        with open(os.path.join(workdir,'readout.temp')) as printout:
            self.printout = printout.read().splitlines()
        os.remove(os.path.join(workdir,'readout.temp'))
        self.output = self.readout()
        return self.output

    def readout(self):
        # Reads off the raw output of the shell command by default, meant to be overloaded
        # Overload this with a function that takes in self.printout and returns something useful to the computation
        print('superclass method called')
        return self.printout

class Rnucleus(Routine):
    def __init__(self,Z,A,neutralMass,I,NDM,NQM):
        """
        Inputs:
        -------
        Z (int): atomic number of nucleus
        A (int): mass number of nucleus
        neutralMass (float): the mass (in amu) of the neutral atom
        I (float): spin of the nucleus given as a decimal
        NDM (float): nuclear dipole moment (nuclear magnetons)
        NQM (float): nuclear quadrupole moment (barns)
        ------
        """
        params = [str(inp) for inp in [Z,A,'n',neutralMass,I,NDM,NQM]]
        super().__init__(name='rnucleus',
                         inputs=[],
                         outputs=['isodata'],
                         params=params)


coredict = {'None':0,'He':1,'Ne':2,'Ar':3,'Kr':4,'Xe':5}
class Rcsfgenerate(Routine):
    def __init__(self,core,csflist,activeset,jlower,jhigher,exc):
        """
        Inputs:
        -------
        core (str): 'None','He','Ne','Ar','Kr','Xe', or 'Rn'
        csflist (list of str): electron configurations
        activeset (list of ints): highest orbital numbers given as a list of quantum numbers in the order s,p,d,f,g,h,etc. So, [4,4,3] means the CSF expansion is truncated at 4s,4p,3d.
        jlower (int): minimum 2*J value of the atom
        jhigher (int): maximum 2*J value of the atom
        exc (int): number of excitations from each config in multireference
        TODO:
        implement default ordering of orbital iteration
        implement multiple different CSF lists with different jlower,jhigher
        ------
        """
        self.header = ['u',str(coredict[core])]
        self.subparams = csflist + ['*']
        self.subparams.append(','.join([str(n)+l for n,l in zip(activeset,['s','p','d','f','g','h','i','l'])]))
        self.subparams.append(f'{jlower},{jhigher}')
        self.subparams.append(str(exc))
        #params.extend(subparams)

        #params.extend(csflist)
        #params.append('*') # end the CSF list
        #params.append(','.join([str(n)+l for n,l in zip(activeset,['s','p','d','f','g','h','i','l'])]))
        #params.append(f'{jlower},{jhigher}')
        #params.append(str(exc))
        #params.append('n') # terminate input
        super().__init__(name='rcsfgenerate',
                         inputs=['clist.ref'],
                         outputs=['rcsf.out','rcsfgenerate.log'],
                         params = self.header + self.subparams + ['n'])
    def __add__(self,other):
        # strip off the terminating 'n' line from self
        # append a 'y' and the subparams from other
        # put the 'n' back
        assert self.header == other.header

        summed = self.__new__(Rcsfgenerate)
        Routine.__init__(summed, name='rcsfgenerate',inputs = ['clist.ref'], outputs = ['rcsf.out','rcsfgenerate.log'],params = self.header + self.subparams + ['y'] + other.subparams + ['n'])
        summed.header = self.header
        summed.subparams = self.subparams + ['y'] + other.subparams
        return summed

    def execute(self,workdir,writeMR=False):
        # we might need to copy to a multiref file at the end
        super().execute(workdir)
        if writeMR:
            copyfile(path.join(workdir,'rcsf.out'),
                    path.join(workdir,'rcsfmr.inp'))
        # then move the .out file to a .inp file for the other functions
        move(path.join(workdir,'rcsf.out'),
             path.join(workdir,'rcsf.inp'))



hamiltoniandict = {'Dirac-Coulomb':1,'Dirac-Coulomb-Breit':2}

class Rcsfinteract(Routine):
    def __init__(self,hamiltonian):
        """
        Inputs:
        -------
        hamiltonian (str): either 'Dirac-Coulomb' or 'Dirac-Coulomb-Breit'
        -----
        """
        #
        params = [str(hamiltoniandict[hamiltonian])]
        super().__init__(name='rcsfinteract',
                         inputs=['rcsfmr.inp','rcsf.inp'],
                         outputs=['rcsf.out'],
                         params=params)

class Rcsfzerofirst(Routine):
    def __init__(self,small_exp,big_exp):
        super().__init__(name='rcsfzerofirst',
                       inputs=[small_exp,big_exp],
                       outputs=['rcsf.out'],
                       params=[small_exp,big_exp])

class Rmixaccumulate(Routine):
    def __init__(self,calcname,useCI,truncate_eps):
        inputs = [f'{calcname}.cm',f'{calcname}.c']
        params = [calcname,booltoyesno(useCI),str(truncate_eps),'y'] # always sort by mixing coeff

        super().__init__(name='rmixaccumulate',
                         inputs=inputs,
                         outputs=['rcsf.out'],
                         params=params)

    def readout(self):
        idx = self.printout.index('         block        ncf') + 1 # start reading iccut on the next line
        iccuttable = pd.DataFrame([line.split() for line in self.printout[idx:]],dtype=float)
        print('called subclass method')
        return [int(val) for val in iccuttable.values[:,1]] #a value for each block

class Rangular(Routine):
    def __init__(self,iccut=None):
        """
        Inputs:
        -------
        None necessary, but you can put in a zero and first order partition ICCUT returned by Rmixaccumulate.execute().
        -------
        """
        if iccut == None:
            params = ['y']
        else:
            params = ['n']
            params.extend([str(e) for e in iccut])
        super().__init__(name='rangular',
                         inputs=['rcsf.inp'],
                         outputs=['rangular.log'],
                         params=params)

methoddict = {'Thomas-Fermi':'2','Screened Hydrogenic':'3'}
class Rwfnestimate(Routine):
    def __init__(self,orbdict=None,fallback='Thomas-Fermi'):
        """
        Inputs:
        -------
        orbdict (dict): a dictionary whose keys are orbital designations (e.g. '5s,5p*' and whose values are relative filepaths from the working directory. If not supplied (not recommended because you will not be able to keep track of the calculation method), we will look for 'rwfn.out' in the working directory.
        fallback (str): Uses the 'Thomas-Fermi' or 'Screened Hydrogenic' method to generate all remaining orbitals
        TODO:
        implement non-default orbital generation parameters
        ------
        """
        if orbdict == None:
            defaultorbs = {'*':'rwfn.out'}
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
        """
        Inputs:
        -------
        asfidx (list of list of ints): list of list of GRASP atomic level serial numbers, block by block
        orbs (list of str): list of orbital designations which are to be varied, e.g. ['5s','5d-','6p*']
        runs (int): maximum number of SCF iterations
        TODO: implement non-default node counting threshold
        implement three grid parameters
        ------
        """
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
        """
        Inputs:
        -------
        name the calculation without a file extension.
        ------
        """
        super().__init__(name=f'rsave {calcname}',
                         inputs = ['rmix.out','rwfn.out','rmcdhf.sum','rmcdhf.log'],
                         outputs= [f'{calcname}.{ext}' for ext in ['c','m','w','sum','log']],
                         params = [])


class Rci(Routine):
    def __init__(self,calcname,includetransverse,modifyfreq,scalefactor,includevacpol,includenms,includesms,estselfenergy,largestn,asfidx):
        """
        Inputs:
        -------
        calcname (str): provide the name of a RMCDHF calculation to use as the basis for an RCI calculation.
        includetransverse (bool)
        modifyfreq (bool)
        scalefactor (str for now): '1.d-6' by default
        includevacpol,includenms,includesms,estselfenergy: (bool)
        largestn (int <= 8)
        asfidx (list of list of ints): see rmcdhf
        ------
        """
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



class Rmixextract(Routine):
    def __init__(self,calcname,useCI,tolerance,sort):
        params = [calcname]
        params.append(booltoyesno(useCI))
        params.append(str(tolerance))
        params.append(booltoyesno(sort))
        print(f'{calcname}.cm')
        super().__init__(name = 'rmixextract',
                    inputs = [f'{calcname}.cm'],
                    outputs= ['rcsf.out'], #really?? shouldn't grasp name it something else
                    params = params)

class JJtoLSJ(Routine):
    def __init__(self,calcname,useCI,unique):
        params = [calcname,booltoyesno(useCI),booltoyesno(unique),'y'] #TODO: implement non-default settings
        inputs = [f'{calcname}.c',f'{calcname}.cm']
        outputs= [f'{calcname}.lsj.lbl'] # and maybe some others
        super().__init__(name = 'jj2lsj',
                    inputs = inputs,
                    outputs= outputs,
                    params = params)

class Rlevels(Routine):
    def __init__(self,calcname):
        params = [f'{calcname}.cm',''] #newline to terminate calcname input
        # TODO: implement multiple calculation results
        super().__init__(name = 'rlevels',
                    inputs = [f'{calcname}.cm'],
                    outputs = [],
                    params = params)
    # TODO: implement multiple calculation results
    def readout(self):
        # reads off the energy levels into a table
        tableheader,start,end = np.where([line[0:2] == '--' for line in self.printout])[0]
        tableList = [line.split() for line in self.printout[start+1:end]]
        print(tableList)
        for line in tableList:
            if len(line) < 8:
                line.extend((8-len(line))*[''])
        df = pd.DataFrame(tableList,columns= ['Number','Position','J','Parity','Energy Total','Levels','Splitting','Configuration'],dtype=float)
        return df



