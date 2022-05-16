"""
author: Ditlev Frickmann
email: reimer.frickmann@gmail.com

This script will create a sigma image for all fits files beneath the specified root directory if they match an input pattern.

"""

from sys import argv
from os.path import isdir
from os import listdir
from re import findall
import numpy as np

from MTLib import files
from MTLib import PATH
from MTLib.fitstools import create_segmentation_map, create_mask_from_segmentation_map


def setup():
    root = argv[1]
    if not isdir(root):
        raise ValueError('Root must be a directory.')
    return root

def main():
    root = setup()
    
    # get fits files in root and below root
    products = ['mask','sigma','segmentation']
    fits_files = [file for file in files.walker(root,extension='.fits') if not np.any([item in file for item in products])]

    print('\n'.join(fits_files))

    for fits_file in fits_files:
        # get name of current source from the fits file
        file_name = files.extract_filename(fits_file)
        
        try:
            print('')
            segmentation_map = create_segmentation_map(fits_file)
            print(f'Segmentation map: {segmentation_map}')
            create_mask_from_segmentation_map(segmentation_map)
        except Exception as e:
            print(f'Failed to create mask for {file_name}: {e}')


if __name__ == '__main__':
    main()