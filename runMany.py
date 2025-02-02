import multiprocessing
from subprocess import call
import os
'''
runName='benchmark_sub2'
dd='/marisdata/petter'
inFiles=[dd+'/pgsOutput'+str(i) for i in [1,2,3]]
L=500
nExp=13
'''
runName='run_sub2'
dd='/marisdata/petter'
inFiles=[dd+'/pgsOutput'+str(i) for i in [1,2,3]]
L=10000
nExp=15

cores=multiprocessing.cpu_count()
wd=os.getcwd();
outFile=dd+'/cgsOutput_'+runName
runFile=wd+'/run_'+runName+'.tmp'
binaryFile=wd+'/bin/cgs_'+runName

call(['gcc',wd+'/cgs.c', '-pthread', '-lfftw3_threads', '-lfftw3', '-lm', '-o',binaryFile])

for Nf in [.1,10]:
    with  open(runFile, 'w') as f:
        f.write(str(len(inFiles))+'\n')
        for inF in inFiles: f.write(inF+'\n')
        f.write(outFile+'\n')
        f.write(str(cores)+'\n')
        f.write(str(Nf)+'\n')
        f.write(str(L)+'\n')
        f.write(str(nExp)+'\n') 
    print('Running Nf='+str(Nf))
    call([binaryFile,runFile])

#call(['rm',binaryFile])
