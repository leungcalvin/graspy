import pdb
import numpy as np
import pandas as pd
import subprocess
from shutil import copyfile,move
import os
from os import path
import warnings

def str_to_float(string):
    if type(string) is not float:
        parts = string.split('D')
        if len(parts) == 1:
            return float(parts[0])
        if len(parts) == 2:
            return float(parts[0]) * 10**(int(parts[1]))

def float_to_str(number):
    return "{:.9e}".format(number).replace('e','D')

def initialize(workdir,clist=None,rm = False):
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
        #print(f'Made path to: {workdir}')
    if rm:
        for old_file in os.scandir(workdir):
            print(f'Removing old file: {old_file}')
            os.unlink(old_file.path)
    if clist is not None:
        with open(os.path.join(workdir,'clist.ref'),'w+') as clistfile:
            # clistfile.write([str(''.join(string+'\n')) for string in clist])
            for c in clist:
                clistfile.write("%s\n" % c)

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
            #print(f'Copying {root}.out -> {root}.inp')
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
        self.hash=hash(self)

    def execute(self,workdir):
        self.workdir = workdir
        inputstr = ''.join(p+'\n' for p in self.params)
        for inp in self.inputs:
            if not checkSmart(path.join(workdir,inp)):
                warnings.warn(f'{self.name}: The input file {inp} is missing')
        process = subprocess.run([f'{self.name} | tee readout.temp'], # capture output to a temp file
                input=inputstr,
                capture_output=False,
                #stdout=subprocess.PIPE,
                shell=True,
                cwd=workdir,
                check=True,
                encoding='utf8')
        missing_files = [outp for outp in self.outputs if not path.exists(path.join(workdir,outp))]
        if len(missing_files) > 0:
            warnings.warn(f'{self.name} failed to produce {missing_files} in {workdir}. This was the command-line input to {self.name}:{self.params}')
        else:
            print(f'{self.name} successfully produced {self.outputs} in {workdir}! This was the command-line input to {self.name}:{self.params}')
        # remove temp file
        with open(os.path.join(workdir,'readout.temp')) as printout:
            self.printout = printout.read().splitlines()
        os.remove(os.path.join(workdir,'readout.temp'))
        self.output = self.readout()
        return self.output

    def readout(self):
        # Reads off the raw output of the shell command by default, meant to be overloaded
        # Overload this with a function that takes in self.printout and returns something useful to the computation
        #print('superclass method called')
        return self.printout

class CSFRoutine(Routine):
    """ A CSFRoutine is any GRASP routine that produces a CSF file as rcsf.out. Since other routines need to read from places like rcsf.inp and rcsfmr.inp, we'll overload the execute() routine to include some file management. We need to make sure that we have a self.write_csf attribute."""
    def execute(self,workdir):
        #print('Calling CSFRoutine execute with directory mgmt!')
        assert hasattr(self,'write_csf'),f'{self.name} has no csf file defined!'
        super().execute(workdir)
        # if we want to keep a log file somewhere, move it over.
        if hasattr(self,'write_log'):
            if self.write_log != 'rcsfgenerate.log':
                move(path.join(workdir,'rcsfgenerate.log'),
                     path.join(workdir,self.write_log))
            # if we need to write to a new file, e.g. a multiref
        if self.write_csf != 'rcsf.out':
            copyfile(path.join(workdir,'rcsf.out'),
                 path.join(workdir,self.write_csf))

        # no matter what, move the .out file to a .inp file for the other functions
        move(path.join(workdir,'rcsf.out'),
             path.join(workdir,'rcsf.inp'))

class MPIRoutine(Routine):
    """
    An MPI Routine implements the MPI version of all the GRASP commands available to us. It makes and sets the MPI_TMP directory at runtime. Do we also need to make a `disks' file?
    """
    def execute_mpi(self,workdir,nproc = 4):
        restore = False
        self.name = f'mpirun -np {nproc} {self.name}_mpi'
        if hasattr(os.environ,'MPI_TMP'):
            restore = True
            old = os.environ['MPI_TMP']
            warnings.warn('Overriding default $MPI_TMP env var.')
        os.environ['MPI_TMP'] = os.path.join(workdir,'mpi_tmp')
        warnings.warn(f'Setting $MPI_TMP : {os.environ["MPI_TMP"]}')
        if not os.path.exists(os.environ['MPI_TMP']):
            warnings.warn(f'Making new directory: {os.environ["MPI_TMP"]}')
            os.makedirs(os.environ['MPI_TMP'])
        #import time
        #time.sleep(5) # wait 5 seconds for directory to be created.
        r = super().execute(workdir)
        if restore:
            os.environ['MPI_TMP'] = old
        return r

class Rnucleus(Routine):
    def __init__(self,Z,A,neutralMass,I,NDM,NQM,rms_radius = None,thickness= None):
        """
        Inputs:
        -------
        Z (int): atomic number of nucleus
        A (int): mass number of nucleus
        neutralMass (float): the mass (in amu) of the neutral atom
        I (float): spin of the nucleus given as a decimal
        NDM (float): nuclear dipole moment (nuclear magnetons)
        NQM (float): nuclear quadrupole moment (barns)
        rms_radius (float): root mean squared radius of nucleus (fm)
        thickness (float): skin thickness of nucleus (fm)

        ------
        """
        if rms_radius is not None and thickness is not None:
            params = [str(inp) for inp in [Z,A,'y',rms_radius,thickness, neutralMass, I, NDM, NQM]]
        else:
            params = [str(inp) for inp in [Z,A,'n',neutralMass,I,NDM,NQM]]
        super().__init__(name='rnucleus',
                         inputs=[],
                         params = params,
                         outputs = ['isodata'])


coredict = {'None':0,'He':1,'Ne':2,'Ar':3,'Kr':4,'Xe':5,'Rn':6}
orderingdict = {'Default':'*','Reverse':'r','Symmetry':'s','User specified':'u'}

class Rcsfgenerate(Routine):
    def __init__(self,core,csflist,active_set,jlower,jhigher,exc,ordering='Default', write_csf = 'rcsf.out'):
        """
        Inputs:
        -------
        core (str): 'None','He','Ne','Ar','Kr','Xe', or 'Rn'
        csflist ((list) of str(s)): electron configurations
        active_set (list of ints): highest orbital numbers given as a list of quantum numbers in the order s,p,d,f,g,h,etc. So, [4,4,3] means the CSF expansion is truncated at 4s,4p,3d.
        jlower (int): minimum 2*J value of the atom
        jhigher (int): maximum 2*J value of the atom
        exc (int): number of excitations from each config in multireference
        TODO:
        implement default ordering of orbital iteration
        implement multiple different CSF lists with different jlower,jhigher
        ------
        """
        if type(csflist) is str:
            csflist = [csflist]
        self.header = [orderingdict[ordering],str(coredict[core])]
        self.subparams = csflist + ['*']
        self.subparams.append(','.join([str(n)+l for n,l in zip(active_set,['s','p','d','f','g','h','i','l'])]))
        self.subparams.append(f'{jlower},{jhigher}')
        self.subparams.append(str(exc))
        self.write_csf = write_csf

        #params.extend(subparams)

        #params.extend(csflist)
        #params.append('*') # end the CSF list
        #params.append(','.join([str(n)+l for n,l in zip(active_set,['s','p','d','f','g','h','i','l'])]))
        #params.append(f'{jlower},{jhigher}')
        #params.append(str(exc))
        #params.append('n') # terminate input
        super().__init__(name='rcsfgenerate',
                         inputs=[],
                         outputs=['rcsf.out','rcsfgenerate.log'],
                         params = self.header + self.subparams + ['n'])
    def __add__(self,other):
        # strip off the terminating 'n' line from self
        # append a 'y' and the subparams from other
        # put the 'n' back
        assert self.header == other.header

        summed = self.__new__(Rcsfgenerate)
        Routine.__init__(summed, name='rcsfgenerate',inputs = [], outputs = ['rcsf.out','rcsfgenerate.log'],params = self.header + self.subparams + ['y'] + other.subparams + ['n'])
        summed.header = self.header
        summed.subparams = self.subparams + ['y'] + other.subparams
        if self.write_csf == other.write_csf:
            summed.write_csf = self.write_csf
        else:
            name,ext=  os.path.splitext(self.write_csf)
            othername,otherext=  os.path.splitext(other.write_csf)
            assert ext == otherext, "Multireferences need to have the same file suffix."

            summed.write_csf = f'{"_".join([name,othername])}{ext}'
        return summed

    def execute(self,workdir):
        super().execute(workdir)
        # if we need to write to a new file, e.g. a multiref
        if self.write_csf != 'rcsf.out':
            copyfile(path.join(workdir,'rcsf.out'),
                 path.join(workdir,self.write_csf))

        # no matter what, move the .out file to a .inp file for the other functions
        move(path.join(workdir,'rcsf.out'),
             path.join(workdir,'rcsf.inp'))

list_to_designation = lambda orblist: ','.join([str(n)+l for n,l in zip(orblist,['s','p','d','f','g','h','i','l'])])

class Rcsfsplit(Routine):
    def __init__(self,calc_name,nsplit,splitorbs,write_csf = 'rcsf.inp'):
        if type(splitorbs[0]) is int: # if splits are specified by maximum principal quantum number value
            splitnames = [str(n) for n in splitorbs]
            splitorbs = [[n]*n for n in splitorbs]
            #print(splitnames,'splitnames')
        else: # generate default splitnames
            splitnames = [str(n) for n in range(nsplit)]
        params = [calc_name,str(nsplit)]
        for splitname,splitorb in zip(splitnames,splitorbs):
            params.append(list_to_designation(splitorb))
            params.append(splitname)
        super().__init__(name=f'rcsfsplit',
                         inputs = [f'{calc_name}.c'],
                         outputs= [f'{calc_name}{ext}.c' for ext in splitnames],
                         params = params)


hamiltoniandict = {'Dirac-Coulomb':1,'Dirac-Coulomb-Breit':2}

class Rcsfinteract(Routine):
    def __init__(self,hamiltonian,write_csf = 'rcsf.out'):
        """
        Inputs:
        -------
        hamiltonian (str): either 'Dirac-Coulomb' or 'Dirac-Coulomb-Breit'
        -----
        """
        #
        self.write_csf = write_csf
        params = [str(hamiltoniandict[hamiltonian])]
        super().__init__(name='rcsfinteract',
                         inputs=['rcsfmr.inp','rcsf.inp'],
                         outputs=['rcsf.out'],
                         params=params)

    def execute(self,workdir):
        super().execute(workdir)
        # if we need to write to a new file, e.g. a multiref
        if self.write_csf != 'rcsf.out':
            copyfile(path.join(workdir,'rcsf.out'),
                 path.join(workdir,self.write_csf))

        # no matter what, move the .out file to a .inp file for the other functions
        move(path.join(workdir,'rcsf.out'),
             path.join(workdir,'rcsf.inp'))

class Rcsfzerofirst(CSFRoutine):
    def __init__(self,small_exp,big_exp,write_csf = 'rcsf.out'):
        self.write_csf = write_csf
        super().__init__(name='rcsfzerofirst',
                       inputs=[small_exp,big_exp],
                       outputs=['rcsf.out'],
                       params=[small_exp,big_exp])

class Rmixaccumulate(Routine):
    def __init__(self,calc_name,use_ci,truncate_eps,write_csf = 'rcsf.out'):
        self.write_csf = write_csf
        inputs = [f'{calc_name}.cm',f'{calc_name}.c']
        params = [calc_name,booltoyesno(use_ci),str(truncate_eps),'y'] # always sort by mixing coeff

        super().__init__(name='rmixaccumulate',
                         inputs=inputs,
                         outputs=['rcsf.out'],
                         params=params)

    def execute(self,workdir):
        super().execute(workdir)
        # if we need to write to a new file, e.g. a multiref
        if self.write_csf != 'rcsf.out':
            copyfile(path.join(workdir,'rcsf.out'),
                 path.join(workdir,self.write_csf))

        # no matter what, move the .out file to a .inp file for the other functions
        move(path.join(workdir,'rcsf.out'),
             path.join(workdir,'rcsf.inp'))

    def readout(self):
        idx = self.printout.index('         block        ncf') + 1 # start reading iccut on the next line
        iccuttable = pd.DataFrame([line.split() for line in self.printout[idx:]],dtype=float)
        #print('called subclass method')
        return [int(val) for val in iccuttable.values[:,1]] #a value for each block

class Rangular(MPIRoutine,Routine):
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
    def __init__(self,orbdict=None,fallback='Thomas-Fermi', grid = None):
        """
        Inputs:
        -------
        orbdict (dict): a dictionary whose keys are orbital designations (e.g. '5s,5p*' and whose values are relative filepaths from the working directory. If not supplied (not recommended because you will not be able to keep track of the calculation method), we will look for 'rwfn.out' in the working directory.
        fallback (str): Uses the 'Thomas-Fermi' or 'Screened Hydrogenic' method to generate all remaining orbitals
        ------
        """
        if grid is None:
            params = ['y']
        else:
            params = ['n','n','','y','n','y',float_to_str(grid['RNT']),float_to_str(grid['H']),str(grid['HP']),str(grid['N'])]
        addtl_files = []
        if orbdict == None:
        # no file input
            defaultorbs = {'*':'rwfn.out'}
        else:
        # get file input, saving wildcard entry for last
            for key in orbdict.keys():
                if key != '*':
                    params.extend(['1',orbdict[key],key])
                    addtl_files.append(orbdict[key])
            if '*' in orbdict.keys():
                params.extend(['1',orbdict['*'],'*'])
                addtl_files.append(orbdict['*'])

        params.extend([methoddict[fallback],'*']) # after looking up everything possible, make sure all orbitals are estimated according to the fallback method
        if grid is not None:
            params.append('n')

        super().__init__(name='rwfnestimate',
                         inputs=['isodata','rcsf.inp']+addtl_files,
                         outputs=['rwfn.inp'],
                         params=params)

WEIGHTINGS = {'Equal':'1','Standard':'5','User [unsupported!]':'9'}
ORTHONORMALIZATIONS = {'Update': '1', 'Self consistency': '2'}
class Rmcdhf(MPIRoutine,Routine):
    def __init__(self,asfidx,orbs,specorbs,runs,weighting_method,grid = None, node_threshold = None, integration_method = None,accuracy = None,orthonormalization_order = 'Update',subruns = None):
        """
        Inputs:
        -------
        asfidx (list of list of ints): list of list of GRASP atomic level serial numbers, block by block
        orbs (list of str): list of orbital designations which are to be varied, e.g. ['5s','5d-','6p*']
        runs (int): maximum number of SCF iterations
        weighting_method ('Equal','Standard','User [not yet supported]'): how to weight different Zeeman levels(?) with a given configuration
        integration_method (1,2,3,4): which non-default integration method to use.
        accuracy: numerical convergence threshold

        TODO: implement non-default node counting threshold

        ------
        """
        # runs is an integer denoting the maximum number of SCF iterations.
        if grid is None and integration_method is None and accuracy is None and orthonormalization_order is not 'Update': # keep default settings.
            params = ['y']
        else: # go deep into the weeds
            params = ['n','n']
            # change the grid
            if grid is not None:
                params.extend(['y','n','y',float_to_str(grid['RNT']),float_to_str(grid['H']),str(grid['HP']),str(grid['N'])])
            else:
                params.extend(['n'])
            # change the accuracy
            if accuracy is None:
                params.extend(['n'])
            else:
                params.extend(['y',float_to_str(accuracy)])
        params.extend( ','.join(str(level) for level in levels) for levels in asfidx)
        params.append(WEIGHTINGS[weighting_method])
        params.append(','.join(orbs))
        params.append(','.join(specorbs))
        params.append(str(runs))
        if integration_method is None:
            params.append('n')
        else:
            params.append('y')
            if node_threshold is None:
               params.append('n')
            else:
                params.extend(['y',float_to_str(node_threshold)])
            if integration_method in [1,2,3,4]:
                params.extend(['y'])
                params.extend(['']*(integration_method - 1))
                params.append('*')
                params.extend(['']*(4 - integration_method))
            else:
                params.extend(['n'])
            # TODO: Implement nondefault options 1) allpositive?, 2) accel parameters for rwfn, 3) accel parameters for evecs, 4) nrefine = 5,
            params.extend(['n','n','n','n'])
            if subruns is int:
                params.extend(['y',str(subruns)])
            else:
                params.extend(['n'])
            # TODO: Currently assumes 4) schmidt orthonormalization = yes
            params.extend(['n'])
        if orthonormalization_order in ORTHONORMALIZATIONS.keys():
            params.append(ORTHONORMALIZATIONS[orthonormalization_order])

        super().__init__(name='rmcdhf',
                         inputs = ['isodata','rcsf.inp','rwfn.inp'],
                         outputs = ['rmix.out','rwfn.out','rmcdhf.sum','rmcdhf.log'],
                         params=params)
    def readout(self):
        try:
            last_index = self.printout.index(" RMCDHF: Execution complete.")
        except ValueError:
            warnings.warn('Yikes. RMCDHF did not converge, or you are using the MPI version for which the convergence checker does not work.')
        start_str = "Subshell    Energy    Method   P0    consistency  Norm-1  factor  JP MTP INV NNP"
        #every time 'Subshell ...' (the headers of the table) is found, a table follows
        start_candidates = [i for i,line in enumerate(self.printout) if line == start_str]
        #print(start_candidates)

        end_str = "Average energy"
        end_candidates = np.array([i for i,line in enumerate(self.printout) if end_str in line])
        #print(end_candidates)
        correct_end_line = end_candidates[end_candidates > start_candidates[-1]][0]

        table = self.printout[start_candidates[-1]+1:correct_end_line-1]
        #prints the contents of the table
        #print(table)
        table = filter(None, table)
        #print(table)
        #removes the line where None is returned
        table_list = [line.split() for line in table]
        #print(table)
        #splits the list of strings of the table into seperate lines
        # print(table_list)
        for line in table_list:
            if len(line) < 11:
                line.extend((11-len(line))*[''])
        df = pd.DataFrame(table_list, columns= ['Subshell', 'Energy','Method','P0','Self-consistency','Norm-1','Damping factor','JP','MTP', "INV", 'NNP'], dtype=float)
        #returns the table, and assigns each column with a header
        # Finally, converts improperly-formatted columns to usable data types.
        for col in ['Energy','P0','Self-consistency', 'Norm-1']:
            df[col] = df[col].apply(str_to_float)
        return df

class Rsave(Routine):
    def __init__(self,calc_name):
        """
        Inputs:
        -------
        name the calculation without a file extension.
        ------
        """
        super().__init__(name=f'rsave {calc_name}',
                         inputs = [],# ['rmix.out','rwfn.out','rmcdhf.sum','rmcdhf.log']; rsave will succeed even if some of these are missing
                         outputs= [f'{calc_name}.{ext}' for ext in ['c','m','w','sum','log']],
                         params = [])

class Rasfsplit(Routine):
    def __init__(self,calc_name,something):
        """
        Inputs:
        -------
        calc_name: name of the calculation without a file extension.
        same (bool): whether csfs are generated from same set of orbitals. Default yes.

        ------
        """
        super().__init__(name=f'rasfsplit',
                         inputs = [f'{calc_name}.c'],
                         outputs= [],
                         params = [booltoyesno(same)])

class Rci(MPIRoutine,Routine):
    def __init__(self,calc_name,include_transverse,modify_freq,scale_factor,include_vacpol,include_nms,include_sms,est_self_energy,largest_n,asfidx):
        """
        Inputs:
        -------
        calc_name (str): provide the name of a RMCDHF calculation to use as the basis for an RCI calculation.
        include_transverse (bool)
        modify_freq (bool)
        scale_factor (str for now): '1.d-6' by default
        include_vacpol,include_nms,include_sms,est_self_energy: (bool)
        largest_n (int <= 8)
        asfidx (list of list of ints): see rmcdhf
        ------
        """
        params = ['y',calc_name]
        params.extend([booltoyesno(param) for param in [include_transverse,modify_freq]])
        if include_transverse:
            params.append(str(scale_factor))
        params.extend([booltoyesno(param) for param in [include_vacpol,include_nms,include_sms,est_self_energy]])
        if est_self_energy:
            params.append(str(largest_n))
        params.extend( ','.join(str(level) for level in levels) for levels in asfidx)
        super().__init__(name = 'rci',
                         inputs = [f'{calc_name}.c',f'{calc_name}.w'],
                         outputs=[f'{calc_name}.cm',f'{calc_name}.csum',f'{calc_name}.clog','rci.res'],
                         params = params)



class Rmixextract(CSFRoutine):
    def __init__(self,calc_name,use_ci,tolerance,sort,write_csf):
        self.write_csf = write_csf
        params = [calc_name]
        params.append(booltoyesno(use_ci))
        params.append(str(tolerance))
        params.append(booltoyesno(sort))
        #print(f'{calc_name}.cm')
        super().__init__(name = 'rmixextract',
                    inputs = [f'{calc_name}.cm'],
                    outputs= ['rcsf.out'], #really?? shouldn't grasp name it something else
                    params = params)

class JJtoLSJ(Routine):
    def __init__(self,calc_name,use_ci,unique):
        params = [calc_name,booltoyesno(use_ci),booltoyesno(unique),'y'] #TODO: implement non-default settings
        inputs = [f'{calc_name}.c']
        if use_ci:
            inputs.append(f'{calc_name}.cm') # only need this if using a CI calculation
        outputs= [f'{calc_name}.lsj.lbl'] # and maybe some others
        super().__init__(name = 'jj2lsj',
                    inputs = inputs,
                    outputs= outputs,
                    params = params)

class Rlevels(Routine):
    def __init__(self,files):
        if type(files) is str:
            files = [files]

        params = files + [''] #newline to terminate calc_name input
        # TODO: implement multiple calculation results
        super().__init__(name = 'rlevels',
                    inputs = list(files),
                    outputs = [],
                    params = params)
    # TODO: implement multiple calculation results
    def readout(self):
        # reads off the energy levels into a table
        tableheader,start,end = np.where([line[0:2] == '--' for line in self.printout])[0]
        tableList = [line.split() for line in self.printout[start+1:end]]
        #print(tableList)
        for line in tableList:
            if len(line) < 8:
                line.extend((8-len(line))*[''])
        df = pd.DataFrame(tableList,columns= ['Number','Position','J','Parity','Energy Total','Levels','Splitting','Configuration'],dtype=float)
        for key in ['Number','Position','J','Energy Total', 'Levels','Splitting']:
            df[key] = pd.to_numeric(df[key])
        return df

class Rhfs(Routine):
    def __init__(self,calc_name,use_ci):
        """
        Inputs:
        calc_name (str): name of calculation without file extension, e.g. 2p_3
        use_ci (bool) : use mixing coefficients from CI calculation?
        """
        params = ['y',calc_name,booltoyesno(use_ci)] #TODO: implement non-default settings
        inputs = ['isodata',f'{calc_name}.w',f'{calc_name}.c',f'{calc_name}.cm']
        outputs= [f'{calc_name}.ch',f'{calc_name}.choffd'] # and maybe some others
        super().__init__(name = 'rhfs',
                    inputs = inputs,
                    outputs= outputs,
                    params = params)

class Rbiotransform(Routine):
    def __init__(self,use_ci,calc_name_initial,calc_name_final,transform_all=True):
        """
        Inputs:
        use_ci (bool) : use mixing coefficients from CI calculation?
        calc_name_initial (str): name of calculation without file extension, e.g. 2s_3
        calc_name_final (str): name of calculation without file extension, e.g. 2p_3. Order of initial/final doesn't matter.
        transform_all (bool): Transform all J symmetries? Default True.
        """
        params = ['y',booltoyesno(use_ci),calc_name_initial,calc_name_final,booltoyesno(transform_all)] #TODO: implement non-default settings
        if calc_name_initial == calc_name_final:
            params.insert(4,'y')
        inputs = ['isodata',f'{calc_name_initial}.c',f'{calc_name_initial}.cm',f'{calc_name_initial}.w',f'{calc_name_final}.c',f'{calc_name_final}.cm',f'{calc_name_final}.w']
        outputs = [f'{calc_name_initial}.cbm',f'{calc_name_initial}.bw',f'{calc_name_initial}.TB',f'{calc_name_final}.cbm',f'{calc_name_final}.bw',f'{calc_name_final}.TB']
        super().__init__(name = 'rbiotransform',
                    inputs = inputs,
                    outputs= outputs,
                    params = params)
class Rtransition(Routine):
    def __init__(self,use_ci,calc_name_initial,calc_name_final,transition_spec):
        """
        Inputs:
        use_ci (bool) : use mixing coefficients from CI calculation?
        calc_name_initial (str): name of calculation without file extension, e.g. 2s_3
        calc_name_final (str): name of calculation without file extension, e.g. 2p_3. Order of initial/final doesn't matter.
        transition_spec (list of str): E.g. ['E1','M2']
        """
        params = ['y',booltoyesno(use_ci),calc_name_initial,calc_name_final, ','.join(transition_spec)] #TODO: implement non-default settings
        inputs = ['isodata',f'{calc_name_final}.w',f'{calc_name_final}.bw',f'{calc_name_final}.cbm',f'{calc_name_initial}.w',f'{calc_name_initial}.bw',f'{calc_name_initial}.cbm',]
        outputs= [f'{calc_name_initial}.{calc_name_final}.ct',f'{calc_name_initial}.{calc_name_final}.-1T']
        super().__init__(name = 'rtransition',
                    inputs = inputs,
                    outputs= outputs,
                    params = params)

class Redf(Routine):
    def __init__(self,use_ci,calc_name):
        """
        Inputs:
        use_ci (bool) : use mixing coefficients from CI calculation?
        calc_name (str): name of calculation without file extension, e.g. 2s_3
        """
        params = ['y',calc_name,booltoyesno(use_ci)] #TODO: implement non-default settings
        inputs = ['isodata',f'{calc_name}.c',f'{calc_name}.w']
        if params[2]:
            inputs.append(f'{calc_name}.cm') # only needed in a CI calculation
            outputs= [f'{calc_name}.ced']
        else:
            outputs= [f'{calc_name}.ed']
        super().__init__(name = 'redf',
                    inputs = inputs,
                    outputs= outputs,
                    params = params)
