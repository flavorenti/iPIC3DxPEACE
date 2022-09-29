import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import glob
import os

## Plot of the scalars output (energy, momentum, charge etc.) of an iPIC3D simulations ##

# species: 0 = electrons 
#	   1 = ions
#	   2 = background electrons
#	   3 = background ions

# path to data
simu_path = sys.argv[1] #'/ccc/scratch/cont005/gen12428/lavorenf/Mercury_SaeInit/PR1/Normal-newbox'
runs = glob.glob( os.path.join(simu_path,'run*') )
nrun = len(runs)

Ncycle = []
TotEnergy = []
TotMomentum = []
EleEnergy = []
MagEnergy = []
KinEnergy = []
KinEnergy = [[],[],[]]
BulEnergy = [[],[],[]]
Charge = [[],[],[]]
Chargei = []

for irun in range(0,nrun):

    if irun==0 :
        sr=1
    else :        
        sr=0

    data = np.loadtxt(simu_path+'/run%i/data/ConservedQuantities.txt'%irun,skiprows=sr)

    ns = (len(data)-6)/3

    Ncycle_1 = data[0] 
    TotEnergy_1 = data[1]
    TotMomentum_1 = data[2]
    EleEnergy_1 = data[3] 
    MagEnergy_1 = data[4]
    KinEnergy_1 = data[5]
    for i in range(0,ns):
    KinEnergy0 = np.concatenate((KinEnergy0,KinEnergy0_1))
    KinEnergy1 = np.concatenate((KinEnergy1,KinEnergy1_1))
    BulEnergy0 = np.concatenate((BulEnergy0,BulEnergy0_1))
    BulEnergy1 = np.concatenate((BulEnergy1,BulEnergy1_1))
    Chargee = np.concatenate((Chargee,Chargee_1))
    Chargei = np.concatenate((Chargei,Chargei_1))

    print(irun)
    print(Ncycle)

# plot evolution conserved quantities
plt.subplot(211)
plt.plot(Ncycle,TotEnergy/TotEnergy[0],label='Total E/E[0]')
plt.plot(Ncycle,TotMomentum/TotMomentum[0],label='Mom/Mom[0]')
plt.plot(Ncycle,Chargee/Chargee[0],label=r'$\rho_e/\rho_e$[0]')
plt.plot(Ncycle,Chargei/Chargei[0],label=r'$\rho_i/\rho_i$[0]')
plt.legend(loc=0,fontsize=11)
plt.yticks(fontsize=13)
plt.xticks([])

plt.subplot(212)
plt.plot(Ncycle,EleEnergy/EleEnergy[0],label=r'$E^2/E^2(0)$')
plt.plot(Ncycle,MagEnergy/MagEnergy[0],label=r'$B^2/B^2(0)$')
plt.plot(Ncycle,KinEnergy/KinEnergy[0],label=r'$P/P(0)$')
plt.plot(Ncycle,KinEnergy0/KinEnergy0[0],label=r'$P_e/P_e(0)$')
plt.plot(Ncycle,KinEnergy1/KinEnergy1[0],label=r'$P_i/P_i(0)$')
plt.plot(Ncycle,BulEnergy0/BulEnergy0[0],label=r'$K_e/K_e(0)$')
plt.plot(Ncycle,BulEnergy1/BulEnergy1[0],label=r'$K_i/K_i(0)$')
plt.legend(loc=0,fontsize=11)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.xlabel('Ncycle',fontsize=14)

plt.tight_layout()
plt.savefig(sys.argv[2]+'/plot_scalars.png',format='png')
plt.close()
