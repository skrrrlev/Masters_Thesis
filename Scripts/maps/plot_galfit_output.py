
from matplotlib import pyplot as plt
import numpy as np
from numpy import exp, log, pi
from re import findall
from astropy import visualization
from astropy.io import fits
from mpl_toolkits.axes_grid1 import make_axes_locatable
from skimage.measure import profile_line
from scipy.special import gamma
from astropy.modeling.models import Sersic2D

name = 's13_F160w'
output_file = 'Output/Maps/s13_F160w/s13_F160w_output.fits'

kappa = lambda n: 1.9992*n - 0.3271
Σe = lambda n, re, ar, mag, zp: ( exp(-(log(10)/2.5)*(mag-zp)) ) / (2*pi*exp(kappa(n))*(kappa(n)**(-2*n))*gamma(2*n)*n*ar*(re**2))
sersic = lambda r, re, n, ar, mag, zp: Σe(n, re, ar, mag, zp) * exp(-kappa(n) * ((r/(re))**(1/n)-1))

def main():

    '''Define the axes and margin sizes'''
    ax_size = 3
    colourbar_size = 0.25
    left_margin = 0.8
    bottom_margin = 0.8
    right_margin = 0.8
    top_margin = 0.2
    inter_margin = 0.1

    '''Get the figure size'''
    fwidth = left_margin + ax_size + inter_margin + ax_size + inter_margin + ax_size + inter_margin + colourbar_size + right_margin
    fheight = top_margin + ax_size + inter_margin + ax_size + bottom_margin

    fig = plt.figure(figsize=(fwidth, fheight),dpi=100)
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
        (bottom_margin+ax_size+inter_margin)/fheight,
        ax_size/fwidth, ax_size/fheight
    ])

    ax1 = fig.add_axes([ 
        (left_margin + ax_size + inter_margin)/fwidth,
        (bottom_margin+ax_size+inter_margin)/fheight,
        ax_size/fwidth, ax_size/fheight
    ])

    ax2 = fig.add_axes([ 
        (left_margin + ax_size + inter_margin + ax_size + inter_margin)/fwidth,
        (bottom_margin+ax_size+inter_margin)/fheight,
        ax_size/fwidth, ax_size/fheight
    ])

    ax_cb = fig.add_axes([ 
        (left_margin + ax_size + inter_margin + ax_size + inter_margin + ax_size + inter_margin)/fwidth,
        (bottom_margin+ax_size+inter_margin)/fheight,
        colourbar_size/fwidth, ax_size/fheight
    ])

    ax3 = fig.add_axes([ 
        (left_margin)/fwidth,
        (bottom_margin)/fheight,
        ax_size/fwidth, ax_size/fheight
    ])

    ax4 = fig.add_axes([ 
        (left_margin + ax_size + inter_margin)/fwidth,
        (bottom_margin)/fheight,
        ax_size/fwidth, ax_size/fheight
    ])

    ax5 = fig.add_axes([ 
        (left_margin + ax_size + inter_margin + ax_size + inter_margin)/fwidth,
        (bottom_margin)/fheight,
        ax_size/fwidth, ax_size/fheight
    ])



    '''Load data'''
    with fits.open(output_file) as hdul:
        image = hdul[1].data
        image_header = hdul[1].header
        model = hdul[2].data
        model_header = hdul[2].header
        residuals = hdul[3].data
    zp = float(image_header['ABMAGZP'])
    
    mask_file, mask_ext = findall(pattern=r'(.+)\[([0-9]+)\]', string=model_header['MASK'])[0]
    with fits.open(mask_file) as hdul:
        mask = hdul[int(mask_ext)].data

    sigma_file, sigma_ext = findall(pattern=r'(.+)\[([0-9]+)\]', string=model_header['SIGMA'])[0]
    with fits.open(sigma_file) as hdul:
        sigma = hdul[int(sigma_ext)].data

    residuals = -100*residuals/image

    '''Get the extent of the map'''
    nrows,ncols = image.shape
    arcsec_per_pixel = float(image_header['CD2_2'])*3600
    xlim = [-ncols/2*arcsec_per_pixel,ncols/2*arcsec_per_pixel]
    ylim = [-nrows/2*arcsec_per_pixel,nrows/2*arcsec_per_pixel]
    extent = (-ncols/2*arcsec_per_pixel,ncols/2*arcsec_per_pixel,-nrows/2*arcsec_per_pixel,nrows/2*arcsec_per_pixel)
    
    '''Plot the image'''
    norm = visualization.simple_norm(image,stretch='linear',percent=95)
    ax0.imshow(image,cmap='bone',norm=norm,origin='lower',extent=extent)
    ax0.axis('equal')
    ax0.set_xlim(xlim)
    ax0.set_ylim(ylim)
    ax0.set_xticklabels([])
    ax0.set_ylabel('Δ (arcsec)',size=12)
    ax0.annotate("Data", xy=(0.01,0.95), xycoords='axes fraction',color='white',size=12,zorder=100)

    '''plot the model'''
    ax1.imshow(model,cmap='bone',norm=norm,origin='lower',extent=extent)
    ax1.set_xlim(xlim)
    ax1.set_ylim(ylim)
    ax1.set_xticklabels([])
    ax1.set_yticklabels([])
    ax1.annotate("Model", xy=(0.01,0.95), xycoords='axes fraction',color='white',size=12,zorder=100)
    
    '''Plot the residuals'''
    residuals_norm = visualization.simple_norm(residuals,stretch='linear',min_cut=-50.0,max_cut=50.0)
    im = ax2.imshow(residuals,cmap='coolwarm',norm=residuals_norm,origin='lower',extent=extent)
    #divider = make_axes_locatable(ax)
    #cax = divider.append_axes('right', size='5%', pad=0.05)
    #fig.colorbar(im, cax=cax, orientation='vertical')
    ax2.set_xlim(xlim)
    ax2.set_ylim(ylim)
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])
    ax2.annotate("Residuals", xy=(0.01,0.95), xycoords='axes fraction',color='white',size=12,zorder=100)

    '''plot the colorbar'''
    cbar = fig.colorbar(im, cax=ax_cb, orientation='vertical')
    cbar.ax.set_ylabel('Percentage', rotation=270, size=12)

    '''Plot the mask'''
    masked_image = np.ma.masked_where(mask, image)
    ax3.imshow(masked_image,cmap='bone',norm=norm,origin='lower',extent=extent)
    ax3.set_xlim(xlim)
    ax3.set_ylim(ylim)
    ax3.set_xlabel('Δ (arcsec)',size=12)
    ax3.set_ylabel('Δ (arcsec)',size=12)
    ax3.annotate("Mask", xy=(0.01,0.95), xycoords='axes fraction',color='white',size=12,zorder=100)

    '''Plot the contours'''
    sigma_sky = 0.0032369466
    levels = np.arange(10,100,20)*sigma_sky
    ax4.imshow(image,cmap='bone',norm=norm,origin='lower',extent=extent)
    ax4.contour(model,levels=levels,extent=extent,cmap='plasma')
    ax4.set_xlim(xlim)
    ax4.set_ylim(ylim)
    ax4.set_yticklabels([])
    ax4.set_xlabel('Δ (arcsec)',size=12)
    ax4.annotate("Contours", xy=(0.01,0.95), xycoords='axes fraction',color='white',size=12,zorder=100)

    '''Plot the 1D model'''
    parameters = ['XC','YC','MAG','RE','N','AR','PA']
    components: list[dict[str,float]] = []
    idx = 1
    while True:
        if 'sky' in model_header[F'COMP_{idx}']:
            break
        component = {par:float(findall(pattern=r'-?[0-9]+\.[0-9]+', string=model_header[f'{idx}_{par}'])[0]) for par in parameters}
        components.append(component)
        idx += 1

    rs = np.linspace(0,extent[1],100)
    for i,c in enumerate(components):
        ys = sersic(rs,c['RE']*arcsec_per_pixel,c['N'],c['AR'],c['MAG'],zp)
        ax5.plot(rs,ys,f'C{i}-',label=f'Component {i+1}')
        ax5.plot(-rs,ys,f'C{i}-')
    ax5.set_yscale('log')
    ax5.legend(frameon=False,prop={'size': 7})
    ax5.yaxis.set_label_position("right")
    ax5.yaxis.tick_right()
    ax5.set_xlim(xlim)
    ax5.set_xlabel('Δ (arcsec)',size=12)  
    ax5.annotate("Components", xy=(0.01,0.95), xycoords='axes fraction',color='black',size=12,zorder=100)  
        
    plt.savefig(f'Figures/galfit_output/{name}.png')

if __name__ == '__main__':
    main()