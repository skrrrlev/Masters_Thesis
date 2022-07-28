# %% Imports and definitions
from matplotlib import pyplot as plt
from matplotlib import patheffects as pe
import numpy as np
from numpy import exp, log, pi
from re import findall
from astropy.visualization import simple_norm
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from scipy.special import gamma
import sep
import math
float_pattern = r'-?[0-9]+\.[0-9]+'
kappa = lambda n: 1.9992*n - 0.3271
Σe = lambda n, re, ar, mag, zp: ( exp(-(log(10)/2.5)*(mag-zp)) ) / (2*pi*exp(kappa(n))*(kappa(n)**(-2*n))*gamma(2*n)*n*ar*(re**2))
sersic = lambda r, re, n, ar, mag, zp: Σe(n, re, ar, mag, zp) * exp(-kappa(n) * ((r/(re))**(1/n)-1))
bone = ['#0B128C','#535473','#2E2F40','#869AA6','#C5D9D7']
light = '#DCE8E8'
dark = '#2E2F40'
stroke = [pe.withStroke(linewidth=2, foreground=dark)]

def slit (field:np.ndarray, theta, xc, yc, ncols, nrows, ps):
    
    theta = (theta+90.0)%360
    
    x = np.arange(0, ncols, 1)
    y = np.arange(0, nrows, 1)
    xx, yy = np.meshgrid(x, y)
    
    xx = xx - xc
    yy = yy - yc
    
    xxrot = xx*math.cos(math.radians(theta))+yy*math.sin(math.radians(theta))
    yyrot = -xx*math.sin(math.radians(theta))+yy*math.cos(math.radians(theta))
    
    
    x2=xxrot[np.abs(np.abs(yyrot)<0.5)]
    
    		
    mean_q=np.zeros(shape=len(x2))
    std16_q=np.zeros(shape=len(x2))
    std84_q=np.zeros(shape=len(x2))	
    	
    path_mask = np.zeros((nrows,ncols),dtype=bool)
    for k in np.arange(0,len(x2)):
        weight = (np.hypot(yyrot, (xxrot-x2[k]))<1)
        w=weight.astype(np.float64)
        w[w<1]=np.nan
        path_mask[w>=1]=True
        
    
    
        mean_q[k]=np.nanmedian(field*w)
        std16_q[k]=np.nanpercentile(field*w, 16.0)
        std84_q[k]=np.nanpercentile(field*w, 84.0)
    
    return x2*ps, mean_q, std16_q, std84_q, path_mask


# %% Variables

name = 's20_F160w'
output_file = f'Output/Maps/{name}/{name}_output.fits'
plot_suffix = '_bd'
norm_percentage = 99
norm_type = 'linear'
residual_levels = np.array([-10,0,10])
forced_pam = None

def main():

    '''Load data'''
    with fits.open(output_file) as hdul:
        image = hdul[1].data
        image_header = hdul[1].header
        wcs = WCS(image_header)
        model = hdul[2].data
        model_header = hdul[2].header
        residuals = hdul[3].data
    
    mask_file, mask_ext = findall(pattern=r'(.+)\[([0-9]+)\]', string=model_header['MASK'])[0]
    with fits.open(mask_file) as hdul:
        mask = hdul[int(mask_ext)].data
    masked_image = np.ma.masked_where(mask, image)

    nrows,ncols = image.shape

    '''Define the axes and margin sizes'''
    hw_ratio = nrows/ncols
    ax_width = 3
    ax_height = ax_width*hw_ratio
    colourbar_size = 0.25
    left_margin = 0.8
    bottom_margin = 0.8
    right_margin = 0.8
    top_margin = 0.2
    inter_margin = 0.1

    '''Get the figure size'''
    fwidth = left_margin + ax_width + inter_margin + ax_width + inter_margin + ax_width + inter_margin + colourbar_size + right_margin
    fheight = top_margin + ax_height + inter_margin + ax_height + bottom_margin

    fig = plt.figure(figsize=(fwidth, fheight),dpi=300)
    fig.patch.set_facecolor('white')
    
    '''
    Add the axes using the add_axes() method, which takes the parameters:
        xmin: Horizontal coordinate of the lower left corner (as fraction).
        ymin: Vertical coordinate of the lower left corner (as fraction).
        dx: Width of the subplot (as fraction).
        dy: Height of the subplot (as fraction).
    '''
    ax0 = fig.add_axes([ 
        (left_margin)/fwidth,
        (bottom_margin+ax_height+inter_margin)/fheight,
        ax_width/fwidth, ax_height/fheight
    ])

    ax1 = fig.add_axes([ 
        (left_margin + ax_width + inter_margin)/fwidth,
        (bottom_margin+ax_height+inter_margin)/fheight,
        ax_width/fwidth, ax_height/fheight
    ])

    ax2 = fig.add_axes([ 
        (left_margin + ax_width + inter_margin + ax_width + inter_margin)/fwidth,
        (bottom_margin+ax_height+inter_margin)/fheight,
        ax_width/fwidth, ax_height/fheight
    ])

    ax_cb = fig.add_axes([ 
        (left_margin + ax_width + inter_margin + ax_width + inter_margin + ax_width + inter_margin)/fwidth,
        (bottom_margin+ax_height+inter_margin)/fheight,
        colourbar_size/fwidth, ax_height/fheight
    ])

    ax3 = fig.add_axes([ 
        (left_margin)/fwidth,
        (bottom_margin)/fheight,
        ax_width/fwidth, ax_height/fheight
    ])

    ax4 = fig.add_axes([ 
        (left_margin + ax_width + inter_margin)/fwidth,
        (bottom_margin)/fheight,
        ax_width/fwidth, ax_height/fheight
    ])

    ax5 = fig.add_axes([ 
        (left_margin + ax_width + inter_margin + ax_width + inter_margin)/fwidth,
        (bottom_margin)/fheight,
        ax_width/fwidth, ax_height/fheight
    ])


    '''Get the extent of the map'''
    arcsec_per_pixel = float(image_header['CD2_2'])*3600
    xlim = [-ncols/2*arcsec_per_pixel,ncols/2*arcsec_per_pixel]
    ylim = [-nrows/2*arcsec_per_pixel,nrows/2*arcsec_per_pixel]
    extent = (-ncols/2*arcsec_per_pixel,ncols/2*arcsec_per_pixel,-nrows/2*arcsec_per_pixel,nrows/2*arcsec_per_pixel)

    '''Define the background'''
    bkg = sep.Background(image.byteswap().newbyteorder(),bw=ncols,bh=nrows)
    print(bkg.globalback,bkg.globalrms)
    globalstd = np.sqrt(bkg.globalrms**2-bkg.globalback**2)
    print(globalstd)
    contour_levels = bkg.globalback+np.array([1,2,4,8,16,32,64])*globalstd*3.0
    residuals_cmap_levels = bkg.globalback+residual_levels*globalstd

    '''Define the north and east vectors and the vector origin in plot units'''
    ncp_world = SkyCoord(ra=0,dec=90,frame='icrs',unit='deg') # North Celestial Pole
    ncp_pixel = wcs.world_to_pixel(ncp_world) # Pixel coordinates of the North Celestial Pole
    vector_base = np.array([ncols*0.8,nrows*0.95]) # Base vector for the north and east vectors in pixels
    arrow_origin = np.array([(vector_base[0]-ncols/2),(vector_base[1]-nrows/2)])*arcsec_per_pixel # Origin of the north and east vectors in plot units
    north_vector = (vector_base - ncp_pixel) # Vector from origin to North Celestial Pole
    north_vector = (arcsec_per_pixel*5)*north_vector/np.linalg.norm(north_vector) # Normalise the north vector to 5 pixels
    rot = np.array([[np.cos(np.pi/2), -np.sin(np.pi/2)], [np.sin(np.pi/2), np.cos(np.pi/2)]])
    east_vector = np.dot(rot,north_vector) # Vector from origin towards East with length of 5 pixels
    
    '''Plot the image'''
    norm = simple_norm(image[np.logical_not(mask)].flatten(),stretch=norm_type,percent=norm_percentage)
    ax0.imshow(image,cmap='bone',norm=norm,origin='lower',extent=extent)
    ax0.set_xlim(xlim)
    ax0.set_ylim(ylim)
    ax0.set_xticklabels([])
    ax0.set_ylabel('Δ (arcsec)',size=12)
    ax0.annotate("Data", xy=(0.01,0.95), xycoords='axes fraction',color=light,path_effects=stroke,size=12,zorder=100)
    
    '''Plot the North and East arrows'''
    ax0.arrow(*arrow_origin, *north_vector, color=light, width=0.025, head_width=0.05, head_length=0.05, zorder=1000)
    ax0.annotate("N", xy=(arrow_origin[0]+north_vector[0],arrow_origin[1]+north_vector[1]-(arcsec_per_pixel*1)), xycoords='data',color=light,path_effects=stroke,size=10,zorder=1000,ha='center',va='top')
    ax0.arrow(*arrow_origin, *east_vector, color=light, width=0.025, head_width=0.05, head_length=0.05, zorder=1000)
    ax0.annotate("E", xy=(arrow_origin[0]+east_vector[0]+(arcsec_per_pixel*1),arrow_origin[1]+east_vector[1]), xycoords='data',color=light,path_effects=stroke,size=10,zorder=1000,ha='left',va='center')

    '''plot the model'''
    ax1.imshow(model,cmap='bone',norm=norm,origin='lower',extent=extent)
    ax1.set_xlim(xlim)
    ax1.set_ylim(ylim)
    ax1.set_xticklabels([])
    ax1.set_yticklabels([])
    ax1.annotate("Model", xy=(0.01,0.95), xycoords='axes fraction',color=light,path_effects=stroke,size=12,zorder=100)
    
    '''Plot the residuals'''
    residuals_norm = simple_norm(residuals_cmap_levels,stretch='linear')
    masked_residuals = np.ma.masked_where(mask, residuals)
    im = ax2.imshow(residuals,cmap='coolwarm',norm=residuals_norm,origin='lower',extent=extent)
    ax2.set_xlim(xlim)
    ax2.set_ylim(ylim)
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])
    ax2.annotate("Residuals", xy=(0.01,0.95), xycoords='axes fraction',color=light,path_effects=stroke,size=12,zorder=100)

    '''plot the colorbar'''
    cbar = fig.colorbar(im, cax=ax_cb, orientation='vertical', ticks=residuals_cmap_levels)
    cbar.ax.set_yticklabels([f'{x:.0f}' for x in (residual_levels)])
    cbar.set_label(r'Residuals [$\sigma_{sky}$]', rotation=270, labelpad=25)
    

    '''Plot the mask'''
    cmap = plt.cm.get_cmap("bone")
    cmap.set_bad('#d0473d')
    ax3.imshow(masked_image,cmap=cmap,norm=norm,origin='lower',extent=extent)
    ax3.set_xlim(xlim)
    ax3.set_ylim(ylim)
    ax3.set_xlabel('Δ (arcsec)',size=12)
    ax3.set_ylabel('Δ (arcsec)',size=12)
    ax3.annotate("Mask", xy=(0.01,0.95), xycoords='axes fraction',color=light,path_effects=stroke,size=12,zorder=100)

    '''Plot the contours'''
    ax4.contour(image, contour_levels,origin='lower',linewidths=1.5,colors='#869AA6',extent=extent,zorder=3)
    CS = ax4.contour(model, contour_levels,origin='lower',linewidths=1.5,colors=dark,extent=extent,zorder=2)
    '''contour labels
    def fmt(x):
        x = x/cont
        s = f"{x:.1f}"
        if s.endswith("0"):
            s = f"{x:.0f}"
        return rf"{s}" if plt.rcParams["text.usetex"] else f"{s}"
    ax4.clabel(CS, v, inline=True, fmt=fmt, fontsize=10,zorder=100)
    '''
    #ax4.imshow(image,cmap='bone',norm=norm,origin='lower',extent=extent)
    #ax4.contour(model,levels=levels,extent=extent,cmap='plasma')
    ax4.set_xlim(xlim)
    ax4.set_ylim(ylim)
    ax4.set_yticklabels([])
    ax4.set_xlabel('Δ (arcsec)',size=12)
    ax4.annotate("Contours", xy=(0.01,0.95), xycoords='axes fraction',color=light,path_effects=stroke,size=12,zorder=100)
    

    '''Plot the 1D profile'''
    xm = float(findall(float_pattern,model_header['1_XC'])[0])-1
    ym = float(findall(float_pattern,model_header['1_YC'])[0])-1
    pam = forced_pam if forced_pam else float(findall(float_pattern,model_header['1_PA'])[0])
    
    r, data_profile_1d, _, _, data_profile_path = slit(image, pam, xm, ym, ncols, nrows, arcsec_per_pixel)
    data_profile_1d[np.isnan(data_profile_1d)] = 0
    ax5.scatter(r,data_profile_1d/globalstd,color='#869AA6',alpha=0.75,edgecolors=dark,zorder=10,label='Data')

    r, model_profile_1d, _, _, model_profile_path = slit(model, pam, xm, ym, ncols, nrows, arcsec_per_pixel)
    model_profile_1d[np.isnan(model_profile_1d)] = 0
    ax5.plot(r,model_profile_1d/globalstd,'-',color=dark,zorder=1,linewidth=1,label='Model')

    ax5.legend(frameon=False,prop={'size': 7})
    ax5.yaxis.set_label_position("right")
    ax5.yaxis.tick_right()
    ax5.set_xlim(*xlim)
    ax5.set_xlabel('Δ (arcsec)',size=12)
    ax5.set_ylabel(r'[$\sigma_{sky}$]',rotation=270,labelpad=25,size=12)
    ax5.annotate("Profile", xy=(0.01,0.95), xycoords='axes fraction',color=light,path_effects=stroke,size=12,zorder=100)
    
    ''' Add path to contour plot '''
    ax4.imshow(data_profile_path,cmap='binary',origin='lower',extent=extent,zorder=1,alpha=0.2)

    ''' Save the figure '''
    plt.savefig(f'Figures/galfit_output/{name}{plot_suffix}.png')
    plt.close()

if __name__ == '__main__':
    main()