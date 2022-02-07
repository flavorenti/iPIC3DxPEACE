import sys
import os
import glob

simu_path = sys.argv[1]
output_path = sys.argv[2]

runs = glob.glob( os.path.join(simu_path,'run*') )[::-1]
nrun = len(runs)
cycle = []
fields_comp = []
fields_comm = []
pcls_comp = []
pcls_comm = []
moments_comp = []
moments_comm = []

for i in range(0,nrun):
    ff = open( glob.glob(os.path.join(runs[i],'RUN*'))[0] ,'r')
    for line in ff.readlines():
        if (len(line.split('======= '))>1):
            cycle.append( int(line.split(' ')[3]) )
        if ( len(line.split('avg_|'))>1 ):
            if line.split(' ')[-1]=='fields\n' :
                fields_comp.append( float(line.split(' ')[-4]) )
                fields_comm.append( float(line.split(' ')[-2]) )
            if line.split(' ')[-1]=='particles\n' :
                pcls_comp.append( float(line.split(' ')[-4]) )
                pcls_comm.append( float(line.split(' ')[-2]) )
            if line.split(' ')[-1]=='moments\n' :
                moments_comp.append( float(line.split(' ')[-4]) )
                moments_comm.append( float(line.split(' ')[-2]) )


dims = [len(cycle),len(fields_comp),len(fields_comm),len(pcls_comp),len(pcls_comm),len(moments_comp),len(moments_comm)]

output = open(output_path+'/computational_time.txt','w')
output.write('# cycle\t fields-computation[s]\t fields-communication[s]\t pcls-comp[s]\t pcls-comm[s]\t moments-comp[s]\t moments-comm[s]\n')

for i in range(0,min(dims)):
    output.write('%i\t %.3f\t %.3f\t %.3f\t %.3f\t %.3f\t %.3f\n'%(cycle[i],fields_comp[i],fields_comm[i],pcls_comp[i],pcls_comm[i],moments_comp[i],moments_comm[i]))

