"""
author: Ditlev Frickmann
email: reimer.frickmann@gmail.com

This script will create a sigma image for all fits files beneath the specified root directory if they match an input pattern.

"""

from sys import argv
from os.path import isdir
from os import listdir
from re import findall
from venv import create

from MTLib import files
from MTLib.fitstools import get_pixel_noise_distribution, create_sigma_image
from MTLib import PATH

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
        cutout_dir = fits_file_path.replace('Data','Output') + file_name + '/'
        # find the cutout images 
        if not isdir(cutout_dir):
            continue

        _, std = get_pixel_noise_distribution(
            fits_file= fits_file,
            plot_file_name= PATH.FIGURES.value + 'noise_distribution/' + file_name + '_noise_distribution.png'
        )

        cutouts = [file for file in listdir(cutout_dir) if file.endswith('.fits') and 'sigma' not in file]
        
        for cutout_file in cutouts:
            create_sigma_image(
                fits_path=cutout_dir + cutout_file,
                std=std,
                output_folder=cutout_dir
            )


if __name__ == '__main__':
    main()