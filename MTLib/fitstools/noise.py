from os import makedirs
from os.path import exists, isdir
from scipy.stats import norm

import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt

from ..files import extract_filename, extract_path

Gaussian = lambda x,mu,sigma: (1 / (sigma * np.sqrt( 2*np.pi )) ) * np.exp(-0.5 * (((x-mu)**2)/(sigma**2)) )

def derive_weight_scale(fits_file: str, noise_pixel_mask_file: str, plot_file_name: str=''):
    '''
    Based on the noise pixels, derive the scale factor of the weights map in the second extension of the fits file
    If the <plot_file_name> is not an empty string, the noise distribution will be histogrammed and saved.
    '''
    with fits.open(fits_file) as hdul:
            sci_data: np.ndarray = hdul[0].data
            wht_data: np.ndarray = hdul[1].data
    with fits.open(noise_pixel_mask_file) as hdul:
        noise_pixel_mask = np.array(hdul[0].data,dtype=bool)
    
    sci_noise_pixels = sci_data[noise_pixel_mask].flatten()
    wht_noise_pixels = wht_data[noise_pixel_mask].flatten()

    distribution = sci_noise_pixels * np.sqrt(wht_noise_pixels)

    mu, std = norm.fit(distribution)

    if plot_file_name:
        plot_dir = extract_path(plot_file_name)
        if not isdir(plot_dir):
            makedirs(plot_dir)
        print('Saving',plot_file_name, sep=' ')

        '''The <data_array> is going to be used for the following analysis. Thus, let's just save a few key properties'''
        maxv = np.max(distribution)
        minv = np.min(distribution)

        '''Create bins for the pixel values'''
        n_bins = int(len(distribution)/2000)
        bins = np.linspace(minv,maxv,n_bins)
        bin_width = (bins[1]-bins[0])
        xbins = bins[:-1] + bin_width/2
        
        counts, _ = np.histogram(distribution,bins) # Entire distribution
        normalization = 1/(np.sum(counts) * bin_width)

        x = np.linspace(minv,maxv,1000)
        y = Gaussian(x,mu,std)

        plt.figure(figsize=(7,4), dpi=150)
        plt.step(xbins,counts*normalization,color='k',alpha=0.4,linewidth=3,label='Real distribution')
        plt.plot(x,y,'k-',label='Fitted Gaussian',alpha=0.8)
        plt.legend(frameon=False)
        plt.grid(alpha=0.25)
        plt.xlabel('Bins of ADU/σ')
        plt.ylabel('Density')
        plt.savefig(plot_file_name,bbox_inches='tight')
    
    return 1/std**2

def create_sigma_image_from_noise_distribution(fits_file: str, noise_pixel_mask_file: str, dead_pixel_mask_file:str, noise_distribution_plot_file: str=''):
    '''
    Based on the weights in the second extension of the fits file, the noise pixel map and the dead pixel mask, create a sigma image.
    If the <noise_distribution_plot_file> is not an empty string, the noise distribution will be histogrammed and saved.
    '''
    output_folder = extract_path(fits_file).replace('Data','Output')
    output_name = extract_filename(fits_file) + '_sigma'
    output_file = output_folder + output_name + '.fits'
    if not exists(output_folder):
        makedirs(output_folder)

    '''First we get the weight scale factor'''
    weight_scale = derive_weight_scale(fits_file,noise_pixel_mask_file,noise_distribution_plot_file)
    
    '''Then we load the dead pixel mask'''
    with fits.open(dead_pixel_mask_file) as dead_hdul:
        dead_mask: np.ndarray = np.array(dead_hdul[1].data == 1, dtype=bool)
        i_dead_mask = np.logical_not(dead_mask)

    '''Now we create the sigma image'''
    # Get the weight image in the fits file
    hdu = fits.open(fits_file)[1]
    hdu.header['EXTVER'] = 'SIGMA'

    # Scale the weights image and then convert to the standard deviation map
    # This is only applied to pixels not in the dead pixel map
    hdu.data[i_dead_mask] = np.sqrt(1/(hdu.data[i_dead_mask]*weight_scale))

    # The weight image could have pixels with value=0, but sigma image are used to devide with, so we set a very large σ-value
    # for pixels that are invalid.
    # They are already filtered out by the dead-mask, but just to be sure!
    hdu.data[dead_mask] = 1E6

    # Now we save the sigma image.
    print(f'Saving σ-image: {output_file}')
    hdu.writeto(output_file, overwrite=True)



if __name__ == '__main__':
    '''
    mean, std = get_pixel_noise_distribution(
        fits_file='Data/Maps/HST/source13/source13_F160w.fits',
        plot_file_name='sigma_img_noise_distribution.png',
        n_bin_factor=100
    )

    print(f'Mean: {mean:.2e}',f'Standard deviation {std:.2e}',sep='\n')

    create_sigma_image(
        fits_path='Data/Maps/HST/source13/cutout_source13_F160w/source13_F160w_galaxy.fits',
        std=std,
        output_folder='.'
    )
    '''
    #derive_weight_scale('Data/Maps/HST/source13/source13_F160w.fits','Output/Maps/HST/source13/source13_F160w_noise.fits')
    create_sigma_image_from_noise_distribution(
        'Data/Maps/HST/source13/source13_F160w.fits',
        'Output/Maps/HST/source13/source13_F160w_noise.fits',
        'Output/Maps/HST/source13/source13_F160w_dead.fits',
        'Figures/noise_distribution/source13_F160w_noise_distribution.png'
    )