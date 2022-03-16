import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys

work_dir = sys.argv[1]

times = np.loadtxt(work_dir+'/texts/computational_time.txt',unpack=True)
cycles = times[0]
cumul = []
cumul_nocomm = []
for ic in range(0,len(cycles)):
    cumul.append(sum((times[1]+times[2]+times[3]+times[4]+times[5]+times[6])[0:ic]))
    cumul_nocomm.append(sum((times[1]+times[3]+times[5])[0:ic]))
plt.figure()
plt.subplot(211)
plt.title('Total times [s]')
plt.plot(cycles,times[1]+times[2],label='fields',color='green')
plt.plot(cycles,times[1],linestyle='--',alpha=0.5,color='green')
plt.plot(cycles,times[3]+times[4],label='particles',color='orange')
plt.plot(cycles,times[3],linestyle='--',alpha=0.5,color='orange')
plt.plot(cycles,times[5]+times[6],label='moments',color='red')
plt.plot(cycles,times[5],linestyle='--',alpha=0.5,color='red')
plt.legend(loc=0,fontsize=12)
plt.xticks([])

plt.subplot(212)
plt.title('Communication times [s]')
plt.plot(cycles,times[2],label='fields',color='green')
plt.plot(cycles,times[4],label='particles',color='orange')
plt.plot(cycles,times[6],label='moments',color='red')
plt.legend(loc=0,fontsize=12)
plt.xlabel('Ncycles',fontsize=14)

plt.savefig(work_dir+'/images/comp_time_percycle-plot.png',format='png')
plt.close()

plt.plot(111)
plt.title('Cumulative time [s]')
plt.plot(cycles,cumul_nocomm,color='blue',linestyle='--',label='without comm')
plt.plot(cycles,cumul,color='blue',label='with comm')
plt.legend(loc=0,fontsize=12)
plt.xlabel('Ncycles',fontsize=14)

plt.tight_layout()
plt.savefig(work_dir+'/images/comp_time_cumulative-plot.png',format='png')
plt.close()

