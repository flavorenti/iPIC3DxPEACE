import numpy as np
import matplotlib.pyplot as plt

## Plot of the scalars output (energy, momentum, charge etc.) of an iPIC3D simulations ##

# species: 0 = electrons 
#	   1 = ions
#	   2 = background electrons
#	   3 = background ions

# path to data
simu_path = '/ccc/scratch/cont005/gen12622/lavorenf/Mercury_SaeInit/PR0'
nrun=1

Ncycle = []
TotEnergy = []
TotMomentum = []
EleEnergy = []
MagEnergy = []
KinEnergy = []
KinEnergy0 = []
KinEnergy1 = []
BulEnergy0 = []
BulEnergy1 = []
Charge = []

for irun in range(0,nrun):

    if (irun==0): 
        sr=1
    else:
        sr=0

    Ncycle_1, TotEnergy_1, TotMomentum_1, EleEnergy_1, MagEnergy_1, KinEnergy_1, KinEnergy0_1, KinEnergy1_1, BulEnergy0_1, BulEnergy1_1, Charge_1 = np.loadtxt(simu_path+'/run_new/data/ConservedQuantities.txt',unpack='True',skiprows=sr)

    Ncycle = np.concatenate((Ncycle,Ncycle_1))
    TotEnergy = np.concatenate((TotEnergy,TotEnergy_1))
    TotMomentum = np.concatenate((TotMomentum,TotMomentum_1))
    EleEnergy = np.concatenate((EleEnergy,EleEnergy_1))
    MagEnergy = np.concatenate((MagEnergy,MagEnergy_1))
    KinEnergy = np.concatenate((KinEnergy,KinEnergy_1))
    KinEnergy0 = np.concatenate((KinEnergy0,KinEnergy0_1))
    KinEnergy1 = np.concatenate((KinEnergy1,KinEnergy1_1))
    BulEnergy0 = np.concatenate((BulEnergy0,BulEnergy0_1))
    BulEnergy1 = np.concatenate((BulEnergy1,BulEnergy1_1))
    Charge = np.concatenate((Charge,Charge_1))

    print(irun)
    print(Ncycle)

# plot evolution conserved quantities
plt.subplot(211)
plt.title(simu_path,fontsize=13)
plt.plot(Ncycle,TotEnergy/TotEnergy[0],label='Total E/E[0]')
plt.plot(Ncycle,TotMomentum/TotMomentum[0],label='Total Mom/Mom[0]')
#plt.plot(Ncycle,(Charge+1.)/(Charge[0]+1.),label='Total (Rho+1)/(Rho[0]+1)')
plt.legend(loc=0,fontsize=14)
plt.yticks(fontsize=14)

plt.subplot(212)
plt.plot(Ncycle,EleEnergy/EleEnergy[0],label=r'$E^2/E^2(0)$')
plt.plot(Ncycle,MagEnergy/MagEnergy[0],label=r'$B^2/B^2(0)$')
plt.plot(Ncycle,KinEnergy/KinEnergy[0],label=r'$P/P(0)$')
plt.plot(Ncycle,KinEnergy0/KinEnergy0[0],label=r'$P_e/P_e(0)$')
plt.plot(Ncycle,KinEnergy1/KinEnergy1[0],label=r'$P_i/P_i(0)$')
plt.plot(Ncycle,BulEnergy0/BulEnergy0[0],label=r'$K_e/K_e(0)$')
plt.plot(Ncycle,BulEnergy1/BulEnergy1[0],label=r'$K_i/K_i(0)$')
plt.legend(loc=0,fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel('Ncycle',fontsize=14)

plt.show()
