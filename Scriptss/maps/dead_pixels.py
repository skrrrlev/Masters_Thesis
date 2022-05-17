"""
author: Ditlev Frickmann
email: reimer.frickmann@gmail.com

This script will create a dead pixel mask from the weight image in the second extension of the fits file.

"""

from sys import argv
from os.path import isdir
from os import makedirs
import numpy as np
from astropy.io import fits

from MTLib import files


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
        # get name of current source from the fits file
        file_name = files.extract_filename(fits_file)
        # define the output folder
        fits_file_path = files.extract_path(fits_file)
        output_dir = fits_file_path.replace('Data','Output')
        if not isdir(output_dir):
            makedirs(output_dir)
        output_file = output_dir + file_name + '_dead.fits'
        
        '''Load the HDU and the WCS of the weights map'''
        wht_hdu = fits.open(fits_file)[1]

        dead_pixel_mask = wht_hdu.data == 0

        data = np.zeros(wht_hdu.data.shape)
        data[dead_pixel_mask] = 1

        wht_hdu.data = data
        wht_hdu.header['EXTVER'] = 'DEADMASK'

        print(f'Saving dead pixel mask: {output_file}')
        wht_hdu.writeto(output_file, overwrite=True)


if __name__ == '__main__':
    main()