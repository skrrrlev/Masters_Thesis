"""
author: Ditlev Frickmann
email: reimer.frickmann@gmail.com

This script will create a sigma image for all fits files beneath the specified root directory if they match an input pattern.

"""

from sys import argv
import numpy as np
from os.path import isfile, isdir
from astropy.io import fits
from scipy.ndimage import binary_dilation

from MTLib import files
from MTLib import PATH
from MTLib.fitstools import create_sigma_image_from_noise_distribution

SOURCE_MASK_PIXEL_DILATION = 7 # pixels

def setup():
    root = argv[1]
    if not isdir(root):
        raise ValueError('Root must be a directory.')
    return root

def main():
    root = setup()
    
    # get fits files in root and below root
    fits_files = files.walker(root,extension='.fits')

    for fits_file in fits_files:

        '''get name and path of the current fits file'''
        file_name = files.extract_filename(fits_file)
        fits_file_path = files.extract_path(fits_file)

        '''Define files looked up or produced by this script'''
        output_dir = fits_file_path.replace('Data','Output')
        dead_mask_file = output_dir + file_name + '_dead.fits'
        segmentation_map = output_dir + file_name + '_segmentation.fits'
        noise_mask_file = output_dir + file_name + '_noise.fits'
        
        '''Make sure that the dead mask file exists and load in the mask'''
        if not isfile(dead_mask_file):
            print(f'Skipped creating sigma image for {file_name}. Could not find dead pixel mask at {dead_mask_file}')
            continue
        
        with fits.open(dead_mask_file) as dead_hdul:
            dead_mask: np.ndarray = np.array(dead_hdul[1].data == 1, dtype=bool)

        '''
        Next, we're making sure that the segmentation map exists. 
        If it does, we create an initial source mask from the map.
        To avoid the few pixels included around the edge of sources, that are affected by the sources, we dilute the mask a little bit.
        '''
        if not isfile(segmentation_map):
            print(f'Skipped creating sigma image for {file_name}. Could not find segmentation map at {segmentation_map}')
            continue

        with fits.open(segmentation_map) as sm_hdul:
            source_mask: np.ndarray = np.array(sm_hdul[1].data != 0, dtype=bool)
        source_mask = binary_dilation(source_mask, iterations=SOURCE_MASK_PIXEL_DILATION)

        '''
        Next, we create the mask of noise pixels, by considering all the pixels
        that are neither included in the dead pixel mask or the source mask.
        Then we save the mask, so that we can check if everything went right.
        '''
        noise_pixel_mask = np.logical_not(np.logical_or(dead_mask, source_mask))

        sci_hdu = fits.open(fits_file)[0]
        sci_hdu.header['EXTVER'] = 'NOISEMASK'

        data = np.zeros(sci_hdu.data.shape)
        data[noise_pixel_mask] = 1
        sci_hdu.data = data
        
        print(f'Saving noise pixel mask: {noise_mask_file}')
        sci_hdu.writeto(noise_mask_file, overwrite=True)

        '''Now that we have the pixel noise mask, we can create the sigma map.'''       
        create_sigma_image_from_noise_distribution(
            fits_file=fits_file,
            noise_pixel_mask_file=noise_mask_file,
            dead_pixel_mask_file=dead_mask_file,
            noise_distribution_plot_file=PATH.FIGURES.value + 'noise_distribution/' + file_name + '_noise_distribution.png'
        )


if __name__ == '__main__':
    main()