import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib import colors as clrs
import scienceplots
# plt.style.use(['science','grid','ieee'])
# plt.style.use(['science','notebook','grid'])
plt.style.use(['science','grid'])
plt.rcParams.update({'font.size': 12,
                    'axes.labelsize': 10,
                    'axes.linewidth': 0.95,
                    'xtick.labelsize' : 10,
                    'ytick.labelsize' : 10,
                    'legend.fontsize' : 10,
                    'legend.fancybox' : False,
                    'legend.edgecolor': 'black'})

# opening default map file
df = pd.read_csv("output/map_dense_13s.txt",sep='\t',header=None)
DF = df.to_numpy()
T = DF[:,0]
t = np.unique(T)
xmin = np.min(t)
xmax = np.max(t)
x_len = len(t)
P = DF[:,1]
p = np.unique(P)
ymin = np.min(p)
ymax = np.max(p)
y_len = len(p)
num_spec = len(DF[0,2:])
z = np.zeros((y_len,x_len))

#fig1,ax1 = plt.subplots(1,1,figsize=(4.92, 4.13), sharex = True, sharey = True)
#fig2,ax2 = plt.subplots(1,1,figsize=(4.92, 4.13), sharex = True, sharey = True)
trace = 1e-8
levels=np.logspace(np.log10(trace),np.log10(1),int(np.abs(np.log10(trace))+1))
lct = ticker.LogLocator(base=10)

for ii in range(num_spec):
    for ip in range(y_len):
        for it in range(x_len):
            z[ip,it] = DF[it+x_len*ip,2+ii]
    print(ii)
    fig2,ax2 = plt.subplots(1,1,figsize=(4.92, 4.13), sharex = True, sharey = True)
    cs1=ax2.contourf(t,p,z,levels=levels,
        vmin=trace,vmax=1,locator=lct,extend='min',
        linestyles = 'solid')
    cl1=ax2.contour(t,p,z,levels=levels)
    fig2.colorbar(cs1,ax=ax2,label=r'$\chi$',ticks=lct,drawedges=True)
    ax2.set_yscale('log')
    plt.show()

'''
trace = 1e-6
nt = 200
fig1,ax1 = plt.subplots(2,2,figsize=(4.92*2, 4.13*2), sharex = True, sharey = True)
fig2,ax2 = plt.subplots(2,2,figsize=(4.92*1, 4.13*1), sharex = True, sharey = True)
xmin, xmax, ymin, ymax = 300, 20000, np.min(pressVecS), np.max(pressVecS)
x = np.linspace(xmin, xmax, nt)
y = np.linspace(ymin, ymax, len(pressVecS))
xl = np.log10(x)
yl = np.log10(y)
xx, yy = np.meshgrid(x,y,indexing='xy')
xxl, yyl = np.meshgrid(xl,yl,indexing='xy')
zzN = 0*xx*yy
zzO = 0*xx*yy
zzNO = 0*xx*yy
zze = 0*xx*yy

for idxp in range(len(pressVecS)):
    print(pathsS[idxp])
    df = pd.read_csv(pathsS[idxp],sep="\t",header=None)
    DF = df.to_numpy()

    for idxt in range(len(DF[:,0])):
        ## N DISS
        #if DF[idxt,2] > trace:
        zzN[idxp,idxt] = DF[idxt,2]
        ## O DISS
        #if DF[idxt,4] > trace:
        zzO[idxp,idxt] = DF[idxt,4]
        ## NO
        #if DF[idxt,5] > trace:
        zzNO[idxp,idxt] = DF[idxt,5]
        ## e
        #if DF[idxt,11] > trace:
        zze[idxp,idxt] = DF[idxt,11]

#levels = np.linspace(trace,1,int(np.abs(np.log10(trace))+1))
levels=np.logspace(np.log10(trace),np.log10(1),int(np.abs(np.log10(trace))+1))
lct = ticker.LogLocator(base=10)

im1=ax1[0][0].imshow(np.log10(zzN), cmap ='Greens',alpha=0.6, extent = [xmin, xmax, ymin, ymax],
        vmin = np.log10(trace), vmax = 0, #np.max(DF[:,4]),
        origin ='lower', aspect='auto', interpolation='none')
fig1.colorbar(im1,ax=ax1[0][0],label=r'$\log \chi_N$')
ax1[0][0].set_title(r'$N$ molar fraction')

cs1=ax2[0][0].contourf(xx,yy,zzN,levels=levels,
        vmin=trace,vmax=1,locator=lct,extend='min',
        linestyles = 'solid')
cl1=ax2[0][0].contour(xx,yy,zzN,levels=levels)
fig2.colorbar(cs1,ax=ax2[0][0],label=r'$\chi_N$',ticks=lct,drawedges=True)
ax2[0][0].set_title(r'$N$ molar fraction')
ax2[0][0].set_yscale('log')


im2=ax1[0][1].imshow(np.log10(zzO), cmap ='Blues',alpha=0.6, extent = [xmin, xmax, ymin, ymax],
        vmin = np.log10(trace), vmax = 0, #np.max(DF[:,4]),
        origin ='lower', aspect='auto', interpolation='none')
ax1[0][1].set_title(r'$O$ molar fraction')
fig1.colorbar(im2,ax=ax1[0][1],label=r'$\log \chi_O$')

cs2=ax2[0][1].contourf(xx,yy,zzO,levels=levels,
        vmin=trace,vmax=1,locator=lct,extend='min')
cl2=ax2[0][1].contour(xx,yy,zzO,levels=levels)
fig2.colorbar(cs2,ax=ax2[0][1],label=r'$\chi_O$',ticks=lct,drawedges=True)
ax2[0][1].set_title(r'$O$ molar fraction')
ax2[0][1].set_yscale('log')

im3=ax1[1][0].imshow(np.log10(zzNO), cmap ='Oranges',alpha=0.6, extent = [xmin, xmax, ymin, ymax],
        vmin = np.log10(trace), vmax = 0, #np.max(DF[:,4]),
        origin ='lower', aspect='auto', interpolation='none')
ax1[1][0].set_title(r'$NO$ molar fraction')
fig1.colorbar(im3,ax=ax1[1][0],label=r'$\log \chi_{NO}$')

cs3=ax2[1][0].contourf(xx,yy,zzNO,levels=levels,
        vmin=trace,vmax=1,locator=lct,extend='min')
cl3=ax2[1][0].contour(xx,yy,zzNO,levels=levels)
fig2.colorbar(cs3,ax=ax2[1][0],label=r'$\chi_{NO}$',ticks=lct,drawedges=True)
ax2[1][0].set_title(r'$NO$ molar fraction')
ax2[1][0].set_yscale('log')

im4=ax1[1][1].imshow(np.log10(zze), cmap ='Purples',alpha=0.6, extent = [xmin, xmax, ymin, ymax],
        vmin = np.log10(trace), vmax = 0, #np.max(DF[:,4]),
        origin ='lower', aspect='auto', interpolation='none')
ax1[1][1].set_title(r'$e^-$ molar fraction')
fig1.colorbar(im4,ax=ax1[1][1],label=r'$\log \chi_{e^-}$')

cs4=ax2[1][1].contourf(xx,yy,zze,levels=levels,
        vmin=trace,vmax=1,locator=lct,extend='min')
cl4=ax2[1][1].contour(xx,yy,zze,levels=levels)
fig2.colorbar(cs4,ax=ax2[1][1],label=r'$\chi_{e^-}$',ticks=lct,drawedges=True)
ax2[1][1].set_title(r'$e^-$ molar fraction')
ax2[1][1].set_yscale('log')
#ax1.set_xlim([xmin,xmax])
#ax1.set_ylim([ymin,ymax])

#ax1[0].set_xlabel(r'Temperature [K]')
#ax1[0].set_ylabel(r'Pressure [Pa]')
fig1.text(0.5, 0.04, r'Temperature [K]', ha='center')
fig1.text(0.04, 0.5, r'Pressure [Pa]', va='center', rotation='vertical')
fig1.savefig("mapim.pdf")
fig2.text(0.5, 0.04, r'Temperature [K]', ha='center')
fig2.text(0.04, 0.5, r'Pressure [Pa]', va='center', rotation='vertical')
fig2.savefig("mapcs.pdf")
'''
'''
#####################################################
fig,ax = plt.subplots(2,1,figsize=(4.92*1, 4.13*2.2))
fig.subplots_adjust(wspace=0.2,hspace=0.2)
# fig.suptitle(r'Starting mixture: $79\% N_2$, $21\% O_2$ - Spieces: $N_2$, $N, O_2$, $O$, $NO$, $NO^+$, $e^-$',fontsize=20)
T_vec = DF[:,0] 
ax[0].plot(T_vec,DF[:,1],lw=2,label=r'$N_2$',color='green',ls='-')
ax[0].plot(T_vec,DF[:,2],lw=2,label=r'$N$',color='green',ls='--')
ax[0].plot(T_vec,DF[:,3],lw=2,label=r'$O_2$',color='blue',ls='-')
ax[0].plot(T_vec,DF[:,4],lw=2,label=r'$O$',color='blue',ls='--')
ax[0].plot(T_vec,DF[:,5],lw=2,label=r'$NO$',color='red',ls='-')
ax[0].plot(T_vec,DF[:,6],lw=2,label=r'$N^{+}$',color='#999900',ls='--')
ax[0].plot(T_vec,DF[:,7],lw=2,label=r'$O^{+}$',color='#00CCCC',ls='--')
ax[0].plot(T_vec,DF[:,8],lw=2,label=r'$NO^{+}$',color='red',ls='--')
ax[0].plot(T_vec,DF[:,9],lw=2,label=r'$N_2^{+}$',color='#999900',ls=':')
ax[0].plot(T_vec,DF[:,10],lw=2,label=r'$O_2^{+}$',color='#00CCCC',ls=':')
ax[0].plot(T_vec,DF[:,11],lw=2,label=r'$e^{-}$',color='#CC6600',ls=':')
ax[0].set_xlabel(r'Temperature, $T$ [K]')
ax[0].set_ylabel(r'Molar Fraction, $\chi_i$')
ax[0].set_xlim([np.min(T_vec),np.max(T_vec)])
ax[0].set_ylim([0,1])
ax[0].legend(loc='upper right')

ax[1].plot(T_vec,np.log10(DF[:,1]),lw=2,label=r'$N_2$',color='green',ls='-')
ax[1].plot(T_vec,np.log10(DF[:,2]),lw=2,label=r'$N$',color='green',ls='--')
ax[1].plot(T_vec,np.log10(DF[:,3]),lw=2,label=r'$O_2$',color='blue',ls='-')
ax[1].plot(T_vec,np.log10(DF[:,4]),lw=2,label=r'$O$',color='blue',ls='--')
ax[1].plot(T_vec,np.log10(DF[:,5]),lw=2,label=r'$NO$',color='red',ls='-')
ax[1].plot(T_vec,np.log10(DF[:,6]),lw=2,label=r'$N^{+}$',color='#999900',ls='--')
ax[1].plot(T_vec,np.log10(DF[:,7]),lw=2,label=r'$O^{+}$',color='#00CCCC',ls='--')
ax[1].plot(T_vec,np.log10(DF[:,8]),lw=2,label=r'$NO^{+}$',color='red',ls='--')
ax[1].plot(T_vec,np.log10(DF[:,9]),lw=2,label=r'$N_2^{+}$',color='#999900',ls=':')
ax[1].plot(T_vec,np.log10(DF[:,10]),lw=2,label=r'$O_2^{+}$',color='#00CCCC',ls=':')
ax[1].plot(T_vec,np.log10(DF[:,11]),lw=2,label=r'$e^{-}$',color='#CC6600',ls=':')
ax[1].set_xlabel(r'Temperature [K]')
ax[1].set_ylabel(r'Molar Fraction, $\log_{10} \chi_i$')
ax[1].set_xlim([np.min(T_vec),np.max(T_vec)])
ax[1].set_ylim([-8,1])
ax[1].legend(loc='upper right')

#fig.show()
fig.savefig('/home/albep/Documents/tesi/report_thermo/gfx/Eq11spec.pdf')
'''
