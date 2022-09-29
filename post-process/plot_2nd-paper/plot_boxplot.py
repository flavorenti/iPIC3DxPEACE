import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('Rates_EII_total_10.txt',delimiter=',')

a_,b_,c_,d_,e_,f_,g_,h_,i_, f_photo = np.loadtxt('Rates_EII.txt',unpack=True,dtype=str)
f_photo = np.array(f_photo,dtype=float)

for i in range(0,10):
	data[i] /= f_photo[i]

ll = ['H','He','O','Na','Mg','Al','Si','K','Ca','Mn']

plt.boxplot(data.T, whis='range', labels=ll)

plt.hlines(1,0,11,linestyle='-',color='red',alpha=0.5)
plt.hlines(0.1,0,11,linestyle='-',color='red',alpha=0.5)

plt.ylabel(r'$\nu_{EII}/\nu_{ph}$',fontsize=22)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.yscale('Log')

plt.show()
