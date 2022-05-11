from re import findall
from os import makedirs
from os import makedirs
from os.path import exists, isdir
from scipy.optimize import minimize

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib import pyplot as plt

from ..files import extract_filename, extract_path

method:str = 'L-BFGS-B'
Gaussian = lambda x,mu,sigma: (1 / (sigma * np.sqrt( 2*np.pi )) ) * np.exp(-0.5 * (((x-mu)**2)/(sigma**2)) )

def get_pixel_noise_distribution(fits_file: str, plot_file_name: str, bright_pixel_cutoff:float=1, fitting_region_std_factor:float=0.4, n_bin_factor:int=100):
    
    '''Load data array from fits file'''
    with fits.open(fits_file) as hdu:
        data = hdu[0].data
    
    '''Flatten the array, sort it from lowest to highest and cutoff <bright_pixel_cutoff> % of the brightest pixels.'''
    sorted_data_array = np.sort(np.array(data).flatten())
    bright_cutoff_index = int(len(sorted_data_array)*(1-bright_pixel_cutoff/100))
    data_array = sorted_data_array[:bright_cutoff_index]

    '''The <data_array> is going to be used for the following analysis. Thus, let's just save a few key properties'''
    maxv = np.max(data_array)
    minv = np.min(data_array)
    std = np.std(data_array)
    median = np.median(data_array)

    '''Create bins for the pixel values'''
    n_bins = int(len(data_array)/n_bin_factor)
    bins = np.linspace(minv,maxv,n_bins)

    '''The right hand side of the distribution is skewed due to pixels that are "contamianted" by astronomical sources.
       The fitting region is thus only gonna encompass the right hand side of the distribution, and a little of the left hand side.
       Also because of the skewed distribution we apply the median instead of the mean, to get an estimate of the peak of the distribution.
    '''
    fitting_region_mask = data_array < ( fitting_region_std_factor * std + median )
    bin_mask = bins < ( fitting_region_std_factor * std + median )

    '''Create the histogram distributions based on the bins, and calculate the integral to get the normalization factor <norm>'''
    full_counts, _ = np.histogram(data_array,bins) # Entire distribution
    norm = np.sum(full_counts)*(bins[1]-bins[0])
    region_bins = bins[bin_mask]
    region_counts, _ = np.histogram(data_array[fitting_region_mask],region_bins) # Fitting distribution

    '''Define the fitting function and the bounds of the fitting parameters of the Gaussian'''
    def tominimize(args):
        y = Gaussian(region_bins[:-1],*args)
        residuals = np.abs(y-region_counts/norm)
        return np.sum(residuals)

    mu_bounds = (median-0.25*std, median+0.25*std)
    std_bounds = (0.1*std, 2*std)

    '''Fit the Gaussian to the distribution and extract the best fit parameters.'''
    fit = minimize(fun=tominimize, x0=(median, std), bounds=(mu_bounds, std_bounds))
    params = fit.x

    if plot_file_name:
        plot_dir = extract_path(plot_file_name)
        if not isdir(plot_dir):
            makedirs(plot_dir)
        print('Saving',plot_file_name, sep=' ')
        plt.figure(figsize=(7,4), dpi=150)
        plt.step(bins[:-1],full_counts/norm,color='k',alpha=0.25,linewidth=3,label='Real distribution')
        plt.step(region_bins[:-1],region_counts/norm,color='k',alpha=0.5,linewidth=3,label='Fitting distribution')
        x = np.linspace(np.min(data_array),np.max(data_array),1000)
        y = Gaussian(x,*params)
        plt.plot(x,y,color='k',label='Fitted gaussian')
        plt.grid(alpha=0.25)
        plt.ylabel('Density')
        plt.xlabel('Pixel value bins')
        plt.legend(frameon=False)
        plt.savefig(plot_file_name,bbox_inches='tight')
    
    return params

def create_sigma_image(fits_path: str, std:float, output_folder: str):
    if not exists(output_folder):
        makedirs(output_folder)
    
    hdu = fits.open(fits_path)[0]
    wcs = WCS(hdu.header)

    # Put the cutout image in the FITS HDU
    hdu.data = std*np.ones(hdu.data.shape,dtype=float)

    # Write the cutout to a new FITS file
    fits_file_name = extract_filename(fits_path)
    cutout_filename = output_folder + '/'*(not output_folder.endswith('/')) + fits_file_name + '_' + 'sigma' + '.fits'
    print(f'Saving {cutout_filename}')
    hdu.writeto(cutout_filename, overwrite=True)


if __name__ == '__main__':
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