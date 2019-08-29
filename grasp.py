import subprocess
from os import path

class Routine(object):
    def __init__(self,name,inputs=[],outputs=[],params=[]):
        self.name = name
        self.inputs=inputs
        self.outputs=outputs
        self.params=params
    def execute(self,workdir):
        inputstr = ''.join(p+'\n' for p in self.params)
        for inp in self.inputs:
            assert path.exists(path.join(workdir,inp)), f'{self.name}: The input file {inp} is missing'
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
            assert(path.exists(path.join(workdir,outp)),f'{self.name}: The output file {outp} is missing')

coredict = {'He':1,'Ne':2,'Xe':5}
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





