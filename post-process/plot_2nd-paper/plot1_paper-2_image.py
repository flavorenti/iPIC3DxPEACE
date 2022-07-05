import fibo as fb
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np



#########################################################################
########### Load the cusps ##############################################
#########################################################################
def load_cusp(input_file_cusp,proj=None):

    # load txt file with values cusp thresold
    phi_c,theta_c,comment = np.loadtxt(input_file_cusp,unpack=True,dtype='str')
    lot_cusp = np.array(phi_c,dtype='float')*24./(2*np.pi)
    lat_cusp = 90.-(np.array(theta_c,dtype='float')*180./np.pi)
    x_cusp=[]
    y_cusp=[]
    for i in range(0,len(phi_c)):
        if lot_cusp[i] not in x_cusp:
            if len(lat_cusp[(phi_c==phi_c[i])*(lat_cusp>0)])>0:
                x_cusp.append(lot_cusp[i])
                y_cusp.append(np.min(lat_cusp[(phi_c==phi_c[i])*(lat_cusp>0)]))
            if len(lat_cusp[(phi_c==phi_c[i])*(lat_cusp<0)])>0:
                x_cusp.append(lot_cusp[i])
                y_cusp.append(np.max(lat_cusp[(phi_c==phi_c[i])*(lat_cusp<0)]))
    x_cusp = np.array(x_cusp)
    y_cusp = np.array(y_cusp)
    x_cuspN = x_cusp[np.where(y_cusp>0)]
    y_cuspN = y_cusp[np.where(y_cusp>0)]
    x_cuspS = x_cusp[np.where(y_cusp<0)]
    y_cuspS = y_cusp[np.where(y_cusp<0)]

    # convert to good units for mollweide proj.
    if proj is not None:
        x_cuspN = (x_cuspN-12.)/12.*np.pi    
        x_cuspS = (x_cuspS-12.)/12.*np.pi    
        y_cuspN = y_cuspN/90.*(np.pi/2.)    
        y_cuspS = y_cuspS/90.*(np.pi/2.) 

    print('CUSP LOADED -->'+input_file_cusp)
    
    return [[x_cuspN,y_cuspN],[x_cuspS,y_cuspS]]



#########################################################################
####### Load Latitude-LT map from prec. file  ###########################
#########################################################################
def load_prec(input_file_prec,proj=None):

    data_prec = np.loadtxt(input_file_prec[0])

    # grid coordinates
    lot = data_prec[0][1:]
    lat = []
    for iy,line in enumerate(data_prec):
        if iy>0:
            lat = np.append(lat,line[0])

    # convert grid to good units for mollweide proj.
    if proj is not None:
        lot = (lot-12.)/12.*np.pi
        lat = lat/90.*(np.pi/2.)

    # create grid    
    Lot,Lat = np.meshgrid(lot,lat)

    # actual data load 1st file
    histo = np.zeros((len(lot),len(lat)))
    for iy,line in enumerate(data_prec):
        if iy>0:
            histo[:,iy-1]=line[1:]

    # if 2nd file exist load it too and MULTIPLY by histo
    if len(input_file_prec)>1 :
        data_prec = np.loadtxt(input_file_prec[1])
        for iy,line in enumerate(data_prec):
            if iy>0:
                if 'Flux' in input_file_prec[1]:
                    histo[:,iy-1] += line[1:]
                elif 'Ener' in input_file_prec[1]:
                    histo[:,iy-1] *= line[1:]

    # print status
    print('GRID LatxLT LOADED -->'+input_file_prec[0])
    if len(input_file_prec)>1 :
        print('MAP LatxLT LOADED WITH %i files'%len(input_file_prec))
    if len(input_file_prec)>1 :
        if 'Flux' in input_file_prec[1]:
            histo /= 2.

    # remove negative values (WHY HAVE THES???)
    histo[histo<0] = 1e-20
    histo[np.isnan(histo)] = 1e-20
    histo[histo==0.] = 1e-20

    # add photoionization rate
    '''
    d_fact = (1/0.47)**2
    i1 = 7
    i2 = 23
    if 'EII-H.txt' in input_file_prec[0]:
        histo[i1:i2,:] += 5e-8*d_fact
    if 'EII-He.txt' in input_file_prec[0]:
        histo[i1:i2,:]+= 7e-8*d_fact
    if 'EII-O.txt' in input_file_prec[0]:
        histo[i1:i2,:] += 2e-7*d_fact
    if 'EII-Mg.txt' in input_file_prec[0]:
        histo[i1:i2,:] += 6e-7*d_fact
    if 'EII-Na.txt' in input_file_prec[0]:
        histo[i1:i2,:] += 7e-6*d_fact
    if 'EII-K.txt' in input_file_prec[0]:
        histo[i1:i2,:] += 2e-5*d_fact
    if 'EII-Mn.txt' in input_file_prec[0]:
        histo[i1:i2,:] += 6e-7*d_fact
    '''

    ff = open('Rates_EII_total.txt','a')
    '''
    ff.write(input_file_prec[0].replace('texts_rates/map_PR1-','').replace('.txt','')+'\t')
    ff.write('%.2e'%np.max(histo)+'\t')
    ff.write('%.2e'%np.min(histo)+'\t')
    ff.write('%.2e'%np.mean(histo)+'\t')
    ff.write('%.2e'%np.std(histo)+'\t')
    ff.write('%.2e'%np.median(histo)+'\t')
    ff.write('%.2e'%np.quantile(histo,0.25)+'\t')
    ff.write('%.2e'%np.quantile(histo,0.50)+'\t')
    ff.write('%.2e'%np.quantile(histo,0.75)+'\t')
    '''
    for ix in range(0,len(histo[:,0])):
        for iy in range(0,len(histo[0,:])):
            ff.write('%.2e,'%histo[ix,iy])
    ff.write('\n')

    return [histo,Lot,Lat]



#########################################################################
########### Put plot in given axis  #####################################
#########################################################################
def plot_data(ax,input_file_prec,input_file_cusp=None,proj=None,xlab=0):

    # open cusp data
    if input_file_cusp is not None:
        cusp = load_cusp(input_file_cusp,proj)
    
    # open prec map data
    prec = load_prec(input_file_prec,proj)

    # define pars plot for two cases
    if len(input_file_prec)==1:
        if 'Flux' in input_file_prec[0]:
            cmm = cm.YlOrBr
            vlim = [1e6,1e9]
            clab = r'$\frac{pcls}{cm^2 s}$'
        else:
            cmm = cm.coolwarm
            vlim = [2e1,4e4]
            clab = r'$eV$'
    else:
        if 'Flux' in input_file_prec[1]:
            cmm = cm.Reds
            vlim = [1e6,1e9]
            clab = r'$\frac{pcls}{cm^2 s}$'
        else:
            cmm = cm.coolwarm
            vlim = [1e8,1e13]
            clab = r'$\frac{eV}{cm^2 s}$'
    if 'XRF' in input_file_prec[0]:
        cmm = cm.inferno
        vlim = [1e-15,1e-12]
        clab = r'$\nu_{XRF}$ [$s^{-1}$]'
    if 'EII' in input_file_prec[0]:
        #cmm = cm.copper
        cmm = cm.jet
        vlim = [5e-8,1e-6]
        clab = r'$\nu_{EII}$ [$s^{-1}$]'

    # plot actual data
    im = ax.pcolormesh(prec[1],prec[2],np.transpose(prec[0]),norm=matplotlib.colors.LogNorm(),vmin=vlim[0],vmax=vlim[1],cmap=cmm)
    
    # plot cusps thresolds
    if input_file_cusp is not None:
        ax.plot(cusp[0][0],cusp[0][1],linestyle='--',marker='',color='white')
        ax.plot(cusp[1][0],cusp[1][1],linestyle='--',marker='',color='white')

    # define extra fancy plot pars
    if proj is not None:
        plt.yticks([-np.pi/2,-np.pi/3,-np.pi/6,0,np.pi/6,np.pi/3],fontsize=14)
        plt.xticks([-np.pi,-np.pi*3/4,-np.pi/2,-np.pi/4.,0,np.pi/4,np.pi/2,np.pi*3/4,np.pi],fontsize=14)
        if not xlab:
            ax.set_xticklabels(['0','3','6','9','12','15','18','21','24'],fontsize=14,color='white')
        else:
            ax.set_xticklabels([])
    else:
        plt.yticks(fontsize=14)
        plt.xticks(fontsize=14)

    plt.grid(alpha=0.6)
    if 'Ca' not in input_file_prec[0]:
        ax.vlines(-np.pi/2.,-np.pi/2,np.pi/2,linestyle='--',color='white')
        ax.vlines(np.pi/2.,-np.pi/2,np.pi/2,linestyle='--',color='white')

    return im,clab,prec




#########################################################################
################# Main of the code  #####################################
#########################################################################
def main(input_file_prec,input_file_cusp=None,proj=None):


    # number of plots
    Nplots = len(input_file_prec)
    if Nplots%2==0 and Nplots<10 :
        Nx = Nplots/2
        Ny = 2
    elif Nplots%2==1 and Nplots<10 :
        Nx = (Nplots+1)/2
        Ny = 2
    elif Nplots>10 :
        Nx = Nplots/3
        Ny = 3

    #extra
    Nx = 5
    Ny = 2

    fig_size = (Nx*6,Ny*4) #8,4 $6,4

    # do the picture
    fig = plt.figure(figsize=fig_size)

    # do the plots
    for ip in range(0,Nplots):
        print('DO PLOTS -->',Ny,Nx,ip+1)
        ax = fig.add_subplot(Ny, Nx, ip+1, projection=proj)
        im,clab,data = plot_data(ax,input_file_prec[ip],input_file_cusp[ip],proj,xlab=ip)
        ll = input_file_prec[ip][0].split('_')[2]
        #ax.set_title(ll,fontsize=12)

        #plt.figure(2)
        #plt.plot(np.mean(data[0],axis=1),label=input_file_prec[ip][0].split('-')[-1].replace('.txt',''))

        if proj is not None:

            # add colorbar1 Flux
            if ip==0 :
                cax1 = fig.add_axes([0.31, 0.11, 0.37, 0.02])
                cb1 = fig.colorbar(im, cax=cax1, orientation='horizontal')
                cb1.set_label(clab,rotation=0,fontsize=23)
                cb1.ax.tick_params(labelsize=22)

            # add colorbar2 Energy
            if ip==4000:
                cax2 = fig.add_axes([0.54, 0.1, 0.37, 0.05])
                cb2 = fig.colorbar(im, cax=cax2, orientation='horizontal')
                cb2.set_label(clab,rotation=0,fontsize=20)
                cb2.ax.tick_params(labelsize=22)
        else:
            if ip>0:
                ax.set_xlabel('Local time [h]',fontsize=16)
            ax.set_ylabel('Latitude [deg]',fontsize=16)
            cb = fig.colorbar(im)
            cb.set_label(clab,rotation=0,fontsize=22)
            cb.ax.tick_params(labelsize=22)
            plt.xlim(0,24)
            plt.ylim(-90,90)

        #plt.xlabel('Local Time',fontsize=18)
        #plt.legend()
        #plt.show()

    # add text
    plt.figtext(0.45,0.93,'Northward IMF',fontsize=33)
    plt.figtext(0.45,0.50,'Southward IMF',fontsize=33)
    plt.figtext(0.15,0.87,'Ca',fontsize=32)
    plt.figtext(0.15,0.45,'Ca',fontsize=32)
    plt.figtext(0.30,0.87,'Si',fontsize=32)
    plt.figtext(0.30,0.45,'Si',fontsize=32)
    plt.figtext(0.45,0.87,'Na',fontsize=32)
    plt.figtext(0.45,0.45,'Na',fontsize=32)
    #plt.figtext(0.50,0.9,'Mg',fontsize=32)
    #plt.figtext(0.61,0.9,'Na',fontsize=32)
    #plt.figtext(0.72,0.9,'K',fontsize=32)
    #plt.figtext(0.91,0.9,'Ca',fontsize=32)
    #plt.figtext(0.60,0.9,'Al',fontsize=32)
    plt.figtext(0.60,0.87,'Fe',fontsize=32)
    plt.figtext(0.60,0.45,'Fe',fontsize=32)
    plt.figtext(0.75,0.87,'O',fontsize=32)
    plt.figtext(0.75,0.45,'O',fontsize=32)
    #plt.figtext(0.75,0.9,'Si',fontsize=32)

    #plt.figtext(0.17,0.45,'5e-8',fontsize=32)
    #plt.figtext(0.28,0.45,'7e-8',fontsize=32)
    #plt.figtext(0.39,0.45,'2e-7',fontsize=32)
    #plt.figtext(0.50,0.45,'6e-7',fontsize=32)
    #plt.figtext(0.61,0.45,'7e-6',fontsize=32)
    #plt.figtext(0.72,0.45,'2e-5',fontsize=32)
    #plt.figtext(0.91,0.45,'7e-5',fontsize=32)
    #plt.figtext(0.6,0.45,'7e-4',fontsize=32)
    #plt.figtext(0.83,0.45,'6e-7',fontsize=32)
    #plt.figtext(0.75,0.45,'2e-5',fontsize=32)

    # save figure and close
    if proj is not None:
        plt.subplots_adjust(top=0.95,bottom=0.1,wspace=0.1, hspace=0.0)
    else:
        plt.subplots_adjust(top=0.87,bottom=0.1,wspace=0.2, hspace=0.4)
    plt.savefig('images/'+save_name+'.png',format='png')
    plt.show()
    plt.close()


# cusps
cuspPR0mBz = 'texts/output_cusps_PR0-minusBz_5000.txt'
cuspPR0pBz = 'texts/output_cusps_PR0-plusBz_6400.txt'
cuspPR1mBz = 'texts/output_cusps_PR1-minusBz_6500.txt'
cuspPR1pBz = 'texts/output_cusps_PR1-plusBz_6000.txt'

# PR0
precPR0mBzFe = 'texts/output_precipitation_PR0-minusBz_Flux-e.txt'
precPR0mBz1Fe = 'texts/output_precipitation_PR0-minusBz_0-100s_Flux-e.txt'
precPR0mBz2Fe = 'texts/output_precipitation_PR0-minusBz_100-1000s_Flux-e.txt'
precPR0mBz3Fe = 'texts/output_precipitation_PR0-minusBz_1000s-inf_Flux-e.txt'
precPR0pBzFe = 'texts/output_precipitation_PR0-plusBz_Flux-e.txt'
precPR0pBz1Fe = 'texts/output_precipitation_PR0-plusBz_0-100s_Flux-e.txt'
precPR0pBz2Fe = 'texts/output_precipitation_PR0-plusBz_100-1000s_Flux-e.txt'
precPR0pBz3Fe = 'texts/output_precipitation_PR0-plusBz_1000s-inf_Flux-e.txt'
#PR1
precPR1mBzFe = 'texts/output_precipitation_PR1-minusBz_Flux-e.txt'
precPR1mBz1Fe = 'texts/output_precipitation_PR1-minusBz_0-100s_Flux-e.txt'
precPR1mBz2Fe = 'texts/output_precipitation_PR1-minusBz_100-1000s_Flux-e.txt'
precPR1mBz3Fe = 'texts/output_precipitation_PR1-minusBz_1000s-inf_Flux-e.txt'
precPR1pBzFe = 'texts/output_precipitation_PR1-plusBz_Flux-e.txt'
precPR1pBz1Fe = 'texts/output_precipitation_PR1-plusBz_0-100s_Flux-e.txt'
precPR1pBz2Fe = 'texts/output_precipitation_PR1-plusBz_100-1000s_Flux-e.txt'
precPR1pBz3Fe = 'texts/output_precipitation_PR1-plusBz_1000s-inf_Flux-e.txt'

#Cross-sections EII 
precPR1pBzH   = 'texts_rates/map_PR1-plusBz_EII-H.txt'
precPR1mBzH   = 'texts_rates/map_PR1-minusBz_EII-H.txt'
precPR1pBzNa  = 'texts_rates/map_PR1-plusBz_EII-Na.txt'
precPR1mBzNa  = 'texts_rates/map_PR1-minusBz_EII-Na.txt'
precPR1pBzK   = 'texts_rates/map_PR1-plusBz_EII-K.txt'
precPR1mBzK   = 'texts_rates/map_PR1-minusBz_EII-K.txt'
precPR1pBzSi  = 'texts_rates/map_PR1-plusBz_EII-Si.txt'
precPR1mBzSi  = 'texts_rates/map_PR1-minusBz_EII-Si.txt'
precPR1pBzMg  = 'texts_rates/map_PR1-plusBz_EII-Mg.txt'
precPR1mBzMg  = 'texts_rates/map_PR1-minusBz_EII-Mg.txt'
precPR1pBzHe  = 'texts_rates/map_PR1-plusBz_EII-He.txt'
precPR1mBzHe  = 'texts_rates/map_PR1-minusBz_EII-He.txt'
precPR1pBzO   = 'texts_rates/map_PR1-plusBz_EII-O.txt'
precPR1mBzO   = 'texts_rates/map_PR1-minusBz_EII-O.txt'
precPR1pBzCa   = 'texts_rates/map_PR1-plusBz_EII-Ca.txt'
precPR1mBzCa  = 'texts_rates/map_PR1-minusBz_EII-Ca.txt'
precPR1pBzAl   = 'texts_rates/map_PR1-plusBz_EII-Al.txt'
precPR1mBzAl  = 'texts_rates/map_PR1-minusBz_EII-Al.txt'
precPR1pBzMn   = 'texts_rates/map_PR1-plusBz_EII-Mn.txt'
precPR1mBzMn  = 'texts_rates/map_PR1-minusBz_EII-Mn.txt'
#Cross-sections XRF 
precPR1pBzCa = 'texts_rates/map_PR1-plusBz_XRF-Ca.txt'
precPR1mBzCa = 'texts_rates/map_PR1-minusBz_XRF-Ca.txt'
precPR1pBzSi = 'texts_rates/map_PR1-plusBz_XRF-Si.txt'
precPR1mBzSi = 'texts_rates/map_PR1-minusBz_XRF-Si.txt'
precPR1pBzFe = 'texts_rates/map_PR1-plusBz_XRF-Fe.txt'
precPR1mBzFe = 'texts_rates/map_PR1-minusBz_XRF-Fe.txt'
precPR1pBzO  = 'texts_rates/map_PR1-plusBz_XRF-O.txt'
precPR1mBzO  = 'texts_rates/map_PR1-minusBz_XRF-O.txt'
precPR1pBzNa = 'texts_rates/map_PR1-plusBz_XRF-Na.txt'
precPR1mBzNa = 'texts_rates/map_PR1-minusBz_XRF-Na.txt'

# map list
ml = []
for simu in [precPR1pBzCa,precPR1pBzSi,precPR1pBzNa,precPR1pBzFe,precPR1pBzO]:
    ml.append([simu])
for simu in [precPR1mBzCa,precPR1mBzSi,precPR1mBzNa,precPR1mBzFe,precPR1mBzO]:
    ml.append([simu])

# cusp list
cl = []
for i in range(5):
    cl.append(cuspPR1pBz) 
for i in range(5):
    cl.append(cuspPR1mBz)

save_name = 'Composite_PR1-good_Flux_0-inf_e6'
save_name = 'Sigma_PR1-ALL_ALL_paper-XRF_1'

main(ml,cl,'mollweide')
