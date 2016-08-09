import multiprocessing
from subprocess import call
import os

runName='2'
dd='/marisdata/petter'
inFiles=[dd+'/pgsOutput'+str(i) for i in [1,2,3]]
L=2000
nExp=16

cores=multiprocessing.cpu_count()
wd=os.getcwd();
outFile=dd+'/cgsOutput_'+runName
runFile=wd+'/run_'+runName

for Nf in [.1,1.,10,.03,.3,3]:
    with  open(runFile, 'w') as f:
        f.write(str(len(inFiles))+'\n')
        for inF in inFiles: f.write(inF+'\n')
        f.write(outFile+'\n')
        f.write(str(cores)+'\n')
        f.write(str(Nf)+'\n')
        f.write(str(L)+'\n')
        f.write(str(nExp)+'\n') 
    print('Running Nf='+str(Nf))
    call([wd+'/cgs',runFile])
