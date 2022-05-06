"""
author: Ditlev Frickmann
email: reimer.frickmann@gmail.com

Based on the stardust fitting results in the <files> list, plot the data points in the SFR-M* space and the ΔMS-redshift space.

"""

from astropy.io import fits
from matplotlib import pyplot as plt
from os import mkdir
from os.path import isdir
from typing import Union
import numpy as np
from matplotlib import cm
from scipy import spatial

from MTLib.plots import annotate_scatter_plot
from MTLib import sfr

root = 'Data/stardust/'
files = [
    'COSMOS/outputs/cosmos/cosmos.fits',
    'ADF/outputs/adf/adf.fits',
    'XMM/outputs/xmm/xmm.fits',
    'HATLAS/outputs/hatlas/hatlas.fits'
]

def theoretical_sfr_sargent(mstar: "Union[float,list[float]]", z: "float")-> "np.ndarray[float]":
    '''From Sargent et al., 2014'''
    mstar = np.array([mstar]).flatten() # solar masses

    N = 9.5E-11 * 10**(-0.21*np.log10(mstar/(5E10))) # yr^-1
    SFR = N * np.exp((2.05*z)/(1+0.16*z**1.54))*mstar # solar masses / yr
    
    return SFR

def theoretical_sfr_schreiber(mstar: "Union[float,list[float]]", z: "float")-> "np.ndarray[float]":
    '''From Sargent et al., 2014'''
    mstar = np.array([mstar]).flatten() # solar masses

    m = np.log10(mstar/1E9)
    r = np.log10(1+z)
    s = np.array([max(0,mm-0.36-2.5*r) for mm in m])
    log_SFR = m-0.5+1.5*r-0.3*(s**2)

    SFR = 10**log_SFR
    
    return SFR

def load_data():
    if not isdir(root+'figures'):
        mkdir(root+'figures')
    ids, z = [], []
    sfr, Δsfr = [], []
    stellar_mass, Δstellar_mass = [], []
    for file in files:
        with fits.open(root+file,'readonly') as f:
            HDU: fits.hdu.table.BinTableHDU = f[1]
            fdata: fits.fitsrec.FITS_rec = HDU.data
        for i, idx in enumerate(fdata.field('id')):
            if idx>30 or fdata.field('mstar')[i] == 0 or fdata.field('sfr')[i] == 0:
                continue
            ids.append(idx)
            z.append(fdata.field('z')[i])
            stellar_mass.append(fdata.field('mstar')[i])
            Δstellar_mass.append(fdata.field('e_mstar')[i])
            sfr.append(fdata.field('sfr')[i])
            Δsfr.append(fdata.field('e_sfr')[i])
    return ids,z,stellar_mass,Δstellar_mass,sfr,Δsfr

def basic_plot(stellar_mass,sfr,z,Δstellar_mass,Δsfr):
    plt.figure(num=1,figsize=(10,5),dpi=100)
    plt.scatter(stellar_mass,sfr,c=z,s=20,vmin=0,vmax=4,cmap='Spectral_r')
    plt.colorbar(label='Redshift')
    plt.errorbar(stellar_mass,sfr,yerr=Δsfr,xerr=Δstellar_mass,fmt='none',capsize=3,capthick=1,ecolor='k',alpha=0.5)
    plt.loglog()
    plt.ylabel(r'SFR [M$_\odot$ yr$^{-1}$]')
    plt.xlabel(r'$M_*$[M$_\odot$]')
    plt.grid(alpha=.2)
    plt.savefig(root+'figures/ms.png')

def overview_plot(stellar_mass,sfr,z,zs_to_plot=[1,2,4],zlim=[1,4],logmstarlim=[9,13.5],dex_scatter=0.2,method:str='sargent'):
    '''Get an overview of the sources in the SFR-M* relation.
    
    - Data:
        1. stellar_mass:    [M_sol]
        2. sfr:             [M_sol/yr]
        3. z:               redshift
    - Theoretical SFR-M* relations 
        4. zs_to_plot:      z of theoretical sfrs to plot
        5. zlim:            limits of the z-colourbar
        6. logmstarlim:     log10 of solarmass limits to plot for the theoretical sfrs
        7. dex_scatter:     scatter to show in [dex] (10^dex_scatter; i.e. 10^(0.3) ~ 2)
        8. method:          "schreiber" or "sargent" 
    '''

    cmap = cm.get_cmap('tab20c')
    mstar = np.logspace(*logmstarlim,50)

    theoretical_sfr_dict = {}
    for z_to_plot in zs_to_plot:
        if method == 'sargent':
            theoretical_sfr_dict[z_to_plot] = theoretical_sfr_sargent(mstar,z_to_plot)
        elif method == 'schreiber':
            theoretical_sfr_dict[z_to_plot] = theoretical_sfr_schreiber(mstar,z_to_plot)
        else:
            raise ValueError('Invalid method; use "sargent" or "schreiber"')

    plt.figure(num=2, figsize=(7,4.2), dpi=150)
    for z_to_plot in zs_to_plot:
        dex = 10**dex_scatter
        z_colour = (z_to_plot-min(zlim))/(max(zlim)-min(zlim))
        colour = cmap(z_colour)
        plt.fill_between(mstar,theoretical_sfr_dict[z_to_plot]*dex,theoretical_sfr_dict[z_to_plot]/dex, alpha=0.5,color=colour)
        plt.plot(mstar, theoretical_sfr_dict[z_to_plot],alpha=1,c=colour,label=r'SFR(M$_\odot$,'+f'z={z_to_plot:.1f})')
        #plt.plot(mstar, theoretical_sfr_dict[z_to_plot]*dex,alpha=0.3,c=colour)
        #plt.plot(mstar, theoretical_sfr_dict[z_to_plot]/dex,alpha=0.3,c=colour)
    plt.scatter(stellar_mass,sfr,c=z,s=20,vmin=min(zlim),vmax=max(zlim),cmap=cmap,edgecolors='k',linewidths=0.5)
    plt.colorbar(label='Redshift')
    plt.loglog()
    plt.legend(frameon=False,prop={'size':6})
    plt.xlim((min(mstar),max(mstar)))
    plt.xlabel(r'Stellar mass $\left[M_\odot\right]$')
    plt.ylabel(r'SFR $\left[M_\odot/yr\right]$')
    plt.grid(alpha=0.2)
    plt.savefig(root+f'figures/ms_overview_{method}.png')
    plt.close()

def group_plot(stellar_mass,sfrs,zs,Δstellar_mass,Δsfrs,groups:"list[tuple]",dex_scatter=0.3,dex_starburst=0.6,logmstarlim=[9,14],method:str='sargent'):
    n_groups = len(groups)
    fig, axes = plt.subplots(ncols=n_groups, nrows=1, sharey=True)
    fig.tight_layout()
    fig.subplots_adjust(top=0.9,bottom=0.2, wspace=0)
    fig.set_size_inches(10,10/n_groups)

    mstar = np.logspace(*logmstarlim,50)
    mstar_ticks = np.logspace(*logmstarlim,4)[1:-1]
    mstar_lim = 10**np.array(logmstarlim)

    sfr_lim = np.array([min(sfrs[sfrs>0])*0.7,max(sfrs)*1.3])
    zlim = [np.floor(min(zs)),np.ceil(max(zs))]
    cmap = cm.get_cmap('tab20c')

    for i in range(n_groups):
        group = groups[i]
        ax = axes[i]

        zmin, zmax = min(group), max(group)

        if zmax == np.inf:
            zmax = max(zs)*1.05
            z_range_str = f'{zmin}' + r' $\leq$ z '
            zcentre = zmin
        else:
            z_range_str = f'{zmin}' + r' $\leq$ z $<$ ' + f'{zmax}'
            zcentre = np.mean(group)


        mask = [True if z >= zmin and z < zmax else False for z in zs]
        zs_in_group = zs[mask]
        mstar_in_group = stellar_mass[mask]
        Δmstar_in_group = Δstellar_mass[mask]
        sfrs_in_group = sfrs[mask]
        Δsfrs_in_group = Δsfrs[mask]
        if method == 'sargent':
            tsfr = theoretical_sfr_sargent(mstar,zcentre)
        elif method == 'schreiber':
            tsfr = theoretical_sfr_schreiber(mstar,zcentre)
        else:
            raise ValueError('Invalid method; use "sargent" or "schreiber"')

        opacity = 0.5
        ax.plot(mstar,tsfr,'k',label=r'SFR(M$_\odot$, z='+f'{zcentre:.2f})',alpha=opacity,zorder=0,linewidth=2)
        ax.plot(mstar,tsfr*10**dex_scatter,'k--',alpha=opacity,zorder=0,linewidth=1)
        ax.plot(mstar,tsfr/10**dex_scatter,'k--',alpha=opacity,zorder=0,linewidth=1)
        ax.plot(mstar,tsfr*10**dex_starburst,'k-.',alpha=opacity,zorder=0,linewidth=1)
        im = ax.scatter(mstar_in_group,sfrs_in_group,c=zs_in_group,s=30,vmin=min(zlim),vmax=max(zlim),cmap=cmap,edgecolors='k',linewidths=1.,zorder=10)
        ax.errorbar(mstar_in_group,sfrs_in_group,yerr=Δsfrs_in_group,xerr=Δmstar_in_group,fmt='none',capsize=2,capthick=0.8,ecolor='k',zorder=5)

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(*mstar_lim)
        ax.set_xticks(mstar_ticks)
        ax.legend(frameon=False,prop={'size': 6})

        if i == 0:
            ax.set_ylabel(r'SFR $\left[M_\odot/yr\right]$')
        ax.set_xlabel(r'Stellar mass $\left[M_\odot\right]$')
        ax.set_title(z_range_str)
        ax.set_ylim(*sfr_lim)
        ax.grid(alpha=0.25)
    
    fig.colorbar(im,ax=axes,label='Redshift')
    plt.savefig(root+f'figures/ms_groups_{method}.png',bbox_inches='tight')

def print_table(ids,mstar,sfrs,z):
    n = len(ids)
    print('ID | Mstar  | SFR    | z')
    print('-'*28)
    print('   | M*     | M*/yr  |  ')
    print('-'*28)
    for i in range(n):
        if sfrs[i]<=0 and mstar[i]<=0:
            print(f'{str(ids[i]):2} | ------ | ------ | {z[i]:.3f}' )
        elif sfrs[i]<=0:
            print(f'{str(ids[i]):2} | {str(np.round(np.log10(mstar[i]),3)):6} | ------ | {z[i]:.3f}' )
        elif mstar[i]<=0:
            print(f'{str(ids[i]):2} | ------ | {str(np.round(np.log10(sfrs[i]),3)):6} | {z[i]:.3f}' )
        else:
            print(f'{str(ids[i]):2} | {str(np.round(np.log10(mstar[i]),3)):6} | {str(np.round(np.log10(sfrs[i]),3)):6} | {z[i]:.3f}' )

    print('~'*28)

def plot_deviation_from_ms(mstar,sfrs,zs,ids,Δmstar,Δsfr,method='schreiber'):
    n = len(mstar)

    Δms = []
    Δms_uncertainty = []
    for i in range(n):
        if method.lower() == 'schreiber':
            tsfr = sfr.schreiber(mstar[i],zs[i])
            Δtsfr = sfr.schreiber_uncertainty(mstar[i],zs[i],Δmstar[i])
        elif method.lower() == 'sargent':
            tsfr = sfr.sargent(mstar[i],zs[i])
            Δtsfr = sfr.sargent_uncertainty(mstar[i],zs[i],Δmstar[i])
        else:
            raise ValueError('Method should be "schreiber" or "sargent".')
        Δms.append(sfrs[i]/tsfr[0])
        Δms_uncertainty.append(np.sqrt( (Δsfr[i]/tsfr[0])**2 + ((sfrs[i]*Δtsfr[0])/tsfr[0]**2)**2 ))

    plt.figure(num=4, figsize=(10,5),dpi=150)
    plt.scatter(zs, Δms, s=20,color='#e0c2a5',edgecolors='k',linewidths=0.5,zorder=10) # #7d9ac9
    Δms = np.array(Δms)
    annotate_scatter_plot(
        strings = [f'{idx}' for idx in ids],
        xs = zs,
        ys = Δms,
        pixel_offset=15,
        min_x_offset=12,
        min_y_offset=0,
        distance_scale_factor=3,
        zorder = 20,
        size = 7
    )
    plt.errorbar(zs,Δms,yerr=Δms_uncertainty,fmt='none',capsize=2,capthick=0.8,ecolor='k',zorder=5)
    
    zlim = [min(zs)*0.8,max(zs)*1.2]
    dex = 0.3
    dex_sb = 0.6
    plt.plot(zlim,[10**dex,10**dex,],'k--',alpha=0.5, zorder=0)
    plt.plot(zlim,[10**-dex,10**-dex],'k--',alpha=0.5, zorder=0)
    plt.plot(zlim,[10**dex_sb,10**dex_sb,],'k-.',alpha=0.5, zorder=0)
    plt.xlim(zlim)
    plt.grid(alpha=0.25)
    plt.yscale('log')
    plt.ylim(0.005,500)
    plt.xlabel('Redshift')
    plt.ylabel(r'$\Delta$ms')
    plt.savefig(root+f'figures/Δms_{method}.png')
    plt.close()


def main():
    ids,z,stellar_mass,Δstellar_mass,sfr,Δsfr = load_data()
    
    print_table(ids,stellar_mass,sfr,z)

    plot_deviation_from_ms(
        mstar=stellar_mass,
        sfrs=sfr,
        zs=z,
        ids=ids,
        Δmstar=Δstellar_mass,
        Δsfr=Δsfr,
        method='schreiber'
    )

    plot_deviation_from_ms(
        mstar=stellar_mass,
        sfrs=sfr,
        zs=z,
        ids=ids,
        Δmstar=Δstellar_mass,
        Δsfr=Δsfr,
        method='sargent'
    )

    return

    overview_plot(
        stellar_mass=stellar_mass,
        sfr=sfr,
        z=z,
        zs_to_plot=[1,2,4],
        zlim=[1,4],
        logmstarlim=[9,14],
        dex_scatter=0.3,
        method='sargent'
    )

    overview_plot(
        stellar_mass=stellar_mass,
        sfr=sfr,
        z=z,
        zs_to_plot=[1,2,4],
        zlim=[1,4],
        logmstarlim=[9,14],
        dex_scatter=0.3,
        method='schreiber'
    )
    
    group_plot(
        stellar_mass=np.array(stellar_mass),
        sfrs=np.array(sfr),
        zs=np.array(z),
        Δstellar_mass = np.array(Δstellar_mass),
        Δsfrs = np.array(Δsfr),
        groups=[(0.5,2),(2,2.5),(2.5,3),(3,4)],
        dex_scatter=0.3,
        dex_starburst=0.6,
        logmstarlim=[10,13],
        method='sargent'
    )

    group_plot(
        stellar_mass=np.array(stellar_mass),
        sfrs=np.array(sfr),
        zs=np.array(z),
        Δstellar_mass = np.array(Δstellar_mass),
        Δsfrs = np.array(Δsfr),
        groups=[(0.5,2),(2,2.5),(2.5,3),(3,4)],
        dex_scatter=0.3,
        dex_starburst=0.6,
        logmstarlim=[10,13],
        method='schreiber'
    )

if __name__ == '__main__':
    main()
    #sfrs()