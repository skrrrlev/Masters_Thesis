
"""
author: Ditlev Frickmann
email: reimer.frickmann@gmail.com

Fit data sources using stardust and plot the results.

"""
import sys
from os.path import isfile
import stardust
import numpy as np
import matplotlib.pyplot as plt
from os import mkdir, listdir
from shutil import rmtree
from os.path import isdir



def help():
    usage = '''
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    usage: python run_stardust.py configfile [testname]
    \tRun from the root directory of this project.\n
    Parameters:
               configfile: the config file.
    (Optional) testname: name of the testrun.
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    '''
    print(usage)
    exit()

def edit_name(ctf: stardust.main.ctf, new_name: str) -> None:
    old_path: str = ctf.config['PATH']
    if isdir(old_path):
        rmtree(old_path[:-1])
    base_dir = old_path.replace('test/','')
    items_in_base_dir = listdir(base_dir)

    i=0
    while True:
        temp_name = new_name + f'_{i}'*bool(i)
        if not temp_name in items_in_base_dir:
            new_name = temp_name
            break
        i += 1
    
    new_path = base_dir + new_name + '/'
    new_fig_path = base_dir + new_name + '/figures/'
    new_sed_path = base_dir + new_name + '/seds/'

    ctf.config['OUTPUT_NAME'] = new_name
    ctf.config['PATH'] = new_path
    mkdir(new_path)
    ctf.config["FIGLOC"] = new_fig_path
    mkdir(new_fig_path)
    ctf.config["sedloc"] = new_sed_path
    mkdir(new_sed_path)

def setup():

    number_of_input_args = len(sys.argv)-1

    if not number_of_input_args:
        help()

    if not isfile(sys.argv[1]) or not sys.argv[1].endswith('.config'):
        help()
    else:
        config_file = sys.argv[1]
        print(f'Found config file: {config_file}')

    # If they're not added, a KeyError is raised. 
    kwargs = {
        'custom_umin_indices': None,
        'custom_qpah_indices': None,
        'custom_gamma': None,
    }

    ctf = stardust.main.ctf(idx=None, config_file=config_file, **kwargs)

    if number_of_input_args > 1:
        new_name = sys.argv[2]
        edit_name(ctf, new_name)
    
    print('\nConfig:')
    for key in ctf.config:
        print(f'\t{key}: {ctf.config[key]}')
    
    return ctf

def cosmos_plots(ctf: stardust.main.ctf):

    ID = ctf.tab['id']

    from MTLib.cosmos import Catalogue
    cc = Catalogue()   
    lp_mass = [10**float(M) for M in cc.get_stellar_masses(ID)]
    ez_mass = [(float(M)) for M in cc.get_property(ID,'ez_mass')] 

    # lp mass plot
    plt.figure(figsize=(10,5),dpi=100)
    minn = min(min(ctf.tab['mstar']),min(lp_mass))*0.95
    maxx = max(max(ctf.tab['mstar']),max(lp_mass))*1.05
    plt.plot([minn,maxx],[minn,maxx],'k--')
    #plt.scatter(mstar,lp_mass,s=10,c='k')
    plt.errorbar(ctf.tab['mstar'],lp_mass,xerr=ctf.tab['e_mstar'],fmt='.',ms=10,mfc='k',capsize=3,capthick=1,ecolor='k',mec='k')
    plt.loglog()
    plt.ylabel(r'lp best $M_*$[M$_\odot$]')
    plt.xlabel(r'Fitted $M_*$[M$_\odot$]')
    plt.xlim([minn,maxx])
    plt.ylim([minn,maxx])
    plt.grid(alpha=.2)
    plt.savefig(ctf.config['FIGLOC']+'lp_mass.png')

    # ez mass plot
    plt.figure(figsize=(10,5),dpi=100)
    minn = min(min(ctf.tab['mstar']),min(ez_mass))*0.95
    maxx = max(max(ctf.tab['mstar']),max(ez_mass))*1.05
    plt.plot([minn,maxx],[minn,maxx],'k--')
    #plt.scatter(mstar,ez_mass,s=10,c='k')
    plt.errorbar(ctf.tab['mstar'],ez_mass,xerr=ctf.tab['e_mstar'],fmt='.',ms=10,mfc='k',capsize=3,capthick=1,ecolor='k',mec='k')
    plt.loglog()
    plt.ylabel(r'ez 500 $M_*$[M$_\odot$]')
    plt.xlabel(r'Fitted $M_*$[M$_\odot$]')
    plt.xlim([minn,maxx])
    plt.ylim([minn,maxx])
    plt.grid(alpha=.2)
    plt.savefig(ctf.config['FIGLOC']+'ez_mass.png')

def plot_sed(ctf: stardust.main.ctf, idx):
    gal_id = ctf.tab['id'][idx]
    irdx = np.int_(ctf.best_ir_idx[idx])

    sed_x = ctf.templ  * (1+ctf.zcat[idx])
    sed_opt = np.sum(ctf.best_coeffs[idx,:-3]*ctf.optical_grid,axis = 1)
    sed_agn = np.sum(ctf.best_coeffs[idx,-3:-1]*ctf.agn_grid,axis = 1)
    sed_ir = ctf.best_coeffs[idx,-1]*ctf.dustnorm*ctf.ir_grid[:,irdx[0],irdx[1],irdx[2]]

    sed_tot = sed_opt + sed_agn + sed_ir

    fig = plt.figure(figsize=(10,5))
    ax = fig.add_subplot(1,1,1)

    ax.fill_between(sed_x,0,sed_opt,color='royalblue',alpha=0.2,label='Stellar')
    ax.fill_between(sed_x,0,sed_agn,color='g',alpha=0.2,label='AGN')
    ax.fill_between(sed_x,0,sed_ir,color='maroon',alpha=0.2,label='Dust')

    sed_radio = ctf.radio_sed(idx,sed_x,alpha=-0.75,use_own_mass=True)
    sed_tot += sed_radio

    ax.plot(sed_x,sed_tot,'k',label='Total',lw=3,alpha=0.6,zorder=10)

    points=((ctf.fnu[idx]/ctf.efnu[idx])>=3)

    ax.errorbar(ctf.wav[idx,points],ctf.fnu[idx,points],yerr=ctf.efnu[idx,points],color='red',fmt='s',capsize=5,capthick=1,ms=12,mfc='white',mew=2,barsabove=True)
    ax.scatter(ctf.wav[idx,~points],(ctf.fnu[idx,~points]+3*ctf.efnu[idx,~points]),marker=r'$\downarrow$',s=300,color='red',zorder=11)

    ax.text(0.02, 0.85,'ID '+str(int(ctf.param[idx,0])), color='k',fontsize=20,transform=ax.transAxes)
    ax.text(0.02, 0.75,r'z = {:.2f}'.format(ctf.zcat[idx]), color='k',fontsize=20,transform=ax.transAxes)

    ax.set_ylabel(r'$f_{\nu}$ [mJy]',fontsize=25)
    ax.set_xlabel(r'$\lambda_{obs}$ $[\mu m]$',fontsize=25)
    ax.legend(fontsize=12)

    ax.set_xlim(.01,10**5.7)
    if not gal_id == 1024221:
        ax.set_ylim(10**-6,10**3)
    else:
        ax.set_ylim(10**-8,10**-2)

    ax.grid(alpha=0.4)
    ax.loglog()

    fig.tight_layout()
    plt.savefig(ctf.config['sedloc']+f'{gal_id}.png')

def main(ctf: stardust.main.ctf):


    n_obj = len(ctf.fnu[:,0]) #Number of objects to fit, leave blank to fit everything
    n_proc = -1 #Number of threads to use, -1 for all available threads, 1 for serial mode
    print('Fitting Catalogue:')
    ctf.fit_catalogue(n_proc=n_proc,n_obj=n_obj)
    print('... done')

    print('Extracting data...')
    mstar = np.log10(ctf.tab['mstar'])
    ID = ctf.tab['id']

    print('Creating figures ...')
    # Create figure
    plt.figure(figsize=(10,5))
    plt.scatter(ctf.tab['mstar'],ctf.tab['sfr_opt'],c=ctf.tab['z'],s=20,vmin=0,vmax=4,cmap='Spectral_r')
    plt.colorbar(label='Redshift')
    plt.errorbar(ctf.tab['mstar'],ctf.tab['sfr_opt'],yerr=ctf.tab['e_sfr_opt'],xerr=ctf.tab['e_mstar'],fmt='none',capsize=3,capthick=1,ecolor='k',alpha=0.5)
    plt.loglog()
    plt.ylabel(r'SFR$_{opt}$ [M$_\odot$ yr$^{-1}$]')
    plt.xlabel(r'$M_*$[M$_\odot$]')
    plt.grid(alpha=.2)
    plt.savefig(ctf.config['FIGLOC']+'ms_optical.png')

    # Create figure
    plt.figure(figsize=(10,5))
    plt.scatter(ctf.tab['mstar'],ctf.tab['sfr'],c=ctf.tab['z'],s=20,vmin=0,vmax=4,cmap='Spectral_r')
    plt.colorbar(label='Redshift')
    plt.errorbar(ctf.tab['mstar'],ctf.tab['sfr'],yerr=ctf.tab['e_sfr'],xerr=ctf.tab['e_mstar'],fmt='none',capsize=3,capthick=1,ecolor='k',alpha=0.5)
    plt.loglog()
    plt.ylabel(r'SFR [M$_\odot$ yr$^{-1}$]')
    plt.xlabel(r'$M_*$[M$_\odot$]')
    plt.grid(alpha=.2)
    plt.savefig(ctf.config['FIGLOC']+'ms.png')

    ctf.save_sed()

    for idx in range(n_obj):
        plot_sed(ctf,idx)
    
    if 'cosmos' in ctf.config['PATH'].lower():
        print('Making COSMOS plots')
        cosmos_plots(ctf)


if __name__ == '__main__':
    ctf = setup()
    main(ctf)