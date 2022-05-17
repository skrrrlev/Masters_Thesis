"""
author: Ditlev Frickmann
email: reimer.frickmann@gmail.com

This script will apply crops to all fits files beneath the specified root directory if they have a corresponding region file.

"""

from sys import argv
from os.path import isdir, isfile

from MTLib import files
from MTLib.fitstools import crop_using_region

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
        # Define the region file
        region_file = fits_file_path + 'regions/' + file_name + '.reg'
        if not isfile(region_file):
            print(f'Skipped creating crops for {file_name}. Could not find region file at {region_file}')
            continue
        
        crop_using_region(
                fits_path=fits_file,
                region_path=region_file,
                output_folder=output_dir,
                output_name=file_name,
                extension=0
            )
        
        dead_pixels = output_dir + file_name + '_dead.fits'
        if isfile(dead_pixels):
            crop_using_region(
                fits_path=dead_pixels,
                region_path=region_file,
                output_folder=output_dir,
                output_name=files.extract_filename(dead_pixels),
                extension=1
            )
        
        segmentation_map = output_dir + file_name + '_segmentation.fits'
        if isfile(segmentation_map):
            crop_using_region(
                fits_path=segmentation_map,
                region_path=region_file,
                output_folder=output_dir,
                output_name=files.extract_filename(segmentation_map),
                extension=1
            )

        sigma_map = output_dir + file_name + '_sigma.fits'
        if isfile(sigma_map):
            crop_using_region(
                fits_path=sigma_map,
                region_path=region_file,
                output_folder=output_dir,
                output_name=files.extract_filename(sigma_map),
                extension=1
            )

if __name__ == '__main__':
    main()