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

    data_prec = np.loadtxt(input_file_prec[0],delimiter=',')

    # grid coordinates
    lot = np.linspace(-180,180,50)
    #lot = np.linspace(-np.pi,np.pi,50)
    lat = np.linspace(-90,90,25)
    #lat = np.linspace(-np.pi/2,np.pi/2,25)

    # create grid    
    Lot,Lat = np.meshgrid(lot,lat)

    # actual data load 1st file
    histo = np.zeros((len(lot),len(lat)))
    for iy,line in enumerate(data_prec):
        if iy>0:
            histo[:,iy-1]=line[:]

    # print status
    print('GRID LatxLT LOADED -->'+input_file_prec[0])
    if len(input_file_prec)>1 :
        print('MAP LatxLT LOADED WITH %i files'%len(input_file_prec))
    if len(input_file_prec)>1 :
        if 'Flux' in input_file_prec[1]:
            histo /= 2.

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
        if 'flux' in input_file_prec[0]:
            cmm = cm.Reds
            vlim = [6e26,1e30]
            clab = r'#pcls'
        else:
            cmm = cm.coolwarm
            vlim = [2e28,4e33]
            clab = r'$eV$'

    # plot actual data
    im = ax.pcolormesh(prec[1],prec[2],np.transpose(prec[0]),norm=matplotlib.colors.LogNorm(),vmin=vlim[0],vmax=vlim[1],cmap=cmm)
    
    # define extra fancy plot pars
    if proj is not None:
        plt.yticks([-np.pi/2,-np.pi/3,-np.pi/6,0,np.pi/6,np.pi/3],fontsize=14)
        plt.xticks([-np.pi*0.98,-np.pi*3/4,-np.pi/2,-np.pi/4.,0,np.pi/4,np.pi/2,np.pi*3/4,np.pi*0.98],fontsize=14)
        if not xlab:
            print('ciao')
            #ax.set_xticklabels(['0','3','6','9','12','15','18','21','24'],fontsize=14)
        else:
            ax.set_xticklabels([])
    else:
        plt.yticks(fontsize=14)
        plt.xticks(fontsize=14)

    plt.grid(alpha=0.6)

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
    Nx = 4
    Ny = 2

    fig_size = (Nx*6,Ny*4) #8,4 $6,4

    # do the picture
    fig = plt.figure(figsize=fig_size)

    # do the plots
    for ip in range(0,Nplots):
        print('DO PLOTS -->',Ny,Nx,ip+1)
        ax = fig.add_subplot(Ny, Nx, ip+1, projection=proj)
        im,clab,data = plot_data(ax,input_file_prec[ip],input_file_cusp[ip],proj,xlab=ip)

        if proj is not None:

            # add colorbar1 Flux
            if ip==0 :
                cax1 = fig.add_axes([0.3, 0.1, 0.45, 0.016])
                cb1 = fig.colorbar(im, cax=cax1, orientation='horizontal')
                cb1.set_label(clab,rotation=0,fontsize=20)
                cb1.ax.tick_params(labelsize=16)

        else:
            if ip>3:
                ax.set_xlabel('Longitude [deg]',fontsize=16)
            if ip==0 or ip==4:
                ax.set_ylabel('Latitude [deg]',fontsize=16)
            cb1 = fig.colorbar(im, orientation='vertical')
            cb1.set_label(clab,rotation=0,fontsize=16)
            cb1.ax.tick_params(labelsize=16)
            plt.xlim(-180,180)
            plt.ylim(-90,90)

    # add text
    plt.figtext(0.45,0.9,'Northward IMF',fontsize=32)
    plt.figtext(0.45,0.45,'Southward IMF',fontsize=32)

    # save figure and close
    if proj is not None:
        plt.subplots_adjust(top=0.95,bottom=0.1,wspace=0.2, hspace=0.0)
    else:
        plt.subplots_adjust(top=0.87,bottom=0.1,wspace=0.3, hspace=0.4)
    plt.savefig('images/'+save_name+'.png',format='png')
    plt.show()
    plt.close()


# PR0
precPR0mBzEn_0_100 = 'energy-zero-hundreds-e-negz.csv'
precPR0pBzEn_0_100 = 'energy-zero-hundreds-e-posz.csv'
precPR0mBzEn_100_1000 = 'energy-hundreds-thousands-e-negz.csv'
precPR0pBzEn_100_1000 = 'energy-hundreds-thousands-e-posz.csv'
precPR0mBzEn_1000_inf = 'energy-thousands-inf-e-negz.csv'
precPR0pBzEn_1000_inf = 'energy-thousands-inf-e-posz.csv'
precPR0mBzEn_0_inf = 'energy-zero-inf-e-negz.csv'
precPR0pBzEn_0_inf = 'energy-zero-inf-e-posz.csv'

# map list
ml = []
for simu in [precPR0pBzEn_0_inf,precPR0pBzEn_0_100,precPR0pBzEn_100_1000,precPR0pBzEn_1000_inf]:
    ml.append([simu.replace('-e-','-i-')])
for simu in [precPR0mBzEn_0_inf,precPR0mBzEn_0_100,precPR0mBzEn_100_1000,precPR0mBzEn_1000_inf]:
    ml.append([simu.replace('-e-','-i-')])

# cusp list
cl = []
for i in range(8):
    cl.append(None) 

save_name = 'Composite_i_Ener_ALL_5'

main(ml,cl,None)
