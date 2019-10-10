import os
import sys

def isCSF(string):
    return '(' in string and ')' in string

def JP(string):
    return string[-2:]

def parseCSFList(filename = 'rcsfmr.inp'):
    with open(filename) as f:
        content = f.readlines()
    content = [x.strip() for x in content]
    while content[0]!='CSF(s):':
        content = content[1:]
    content = content[1:]
    # strips off list of orbitals which vary from CSF to CSF list
    print(content[0:3])
    MR = dict()

    lineidx = 0
    while lineidx != len(content):
        line = content[lineidx]
        if isCSF(line):
            JPblock = JP(content[lineidx+2])
            if JPblock in MR.keys():
                MR[JPblock].append((content[lineidx],content[lineidx+1],content[lineidx+2]))
            else:
                MR[JPblock] = [(content[lineidx],content[lineidx+1],content[lineidx+2])]
            lineidx = lineidx+3
        else:
            lineidx = lineidx+1
    return MR

def indexCSFList(queryCSFs,databaseCSFs):
    JPblocks = list(queryCSFs.keys())
    JPblocks.sort()
    asfidxes = []
    for block in JPblocks:
        currentblockidx = []
        for CSF in queryCSFs[block]:
            try:
                currentblockidx.append(databaseCSFs[block].index(CSF)+1) # grasp indexes from one
            except ValueError:
                print(f'The CSF {CSF} of parity {block} was not found in the database, continuing...')

        asfidxes.append(currentblockidx)
        print('ASF Block: '+block)
        print(','.join([str(sn) for sn in currentblockidx]))
    return asfidxes

def extractMRindices(workdir,mrfilename='rcsfmr.inp',expfilename='rcsf.inp'):
    query = parseCSFList(os.path.join(workdir,mrFilename))
    print(query)
    database = parseCSFList(os.path.join(workdir,expFilename))
    asfidxes = indexCSFList(query,database)


if __name__=='__main__':
    mrFilename = sys.argv[1]
    expFilename = sys.argv[2]

    query = parseCSFList(mrFilename)
    database = parseCSFList(expFilename)
    asfidxes = indexCSFList(query,database)
    print(asfidxes)



