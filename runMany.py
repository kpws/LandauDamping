import multiprocessing
from subprocess import call
import os

runName='5'
dd='/marisdata/petter'
inFiles=[dd+'/pgsOutput'+str(i) for i in [1,2,3]]
L=1000
nExp=15

cores=multiprocessing.cpu_count()
wd=os.getcwd();
outFile=dd+'/cgsOutput_'+runName
runFile=wd+'/run_'+runName+'.tmp'

for Nf in [0.001]:
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
