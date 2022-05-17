"""
author: Ditlev Frickmann
email: reimer.frickmann@gmail.com

This script will create a sigma image for all fits files beneath the specified root directory if they match an input pattern.

"""

from sys import argv
from os.path import isdir
from os import listdir
from os import listdir

from MTLib import files
from MTLib.fitstools import create_mask_from_segmentation_map


def setup():
    root = argv[1]
    if not isdir(root):
        raise ValueError('Root must be a directory.')
    return root

def main():
    root = setup()
    if 'Output/' in root:
        raise ValueError('Specify root in data directory.')
    
    # get fits files in root and below root
    fits_files = [file for file in files.walker(root,extension='.fits')]

    for fits_file in fits_files:
        # get name of current source from the fits file
        file_name = files.extract_filename(fits_file)
        output_dir = files.extract_path(fits_file).replace('Data','Output')

        cutout_directories = [output_dir+item+'/' for item in listdir(output_dir) if isdir(output_dir+item+'/')]

        for subdir in cutout_directories:
        
            try:
                segmentation_map = subdir + file_name + '_segmentation.fits'
                create_mask_from_segmentation_map(segmentation_map)
            except Exception as e:
                print(f'Failed to create mask for {segmentation_map}: {e}')

if __name__ == '__main__':
    main()