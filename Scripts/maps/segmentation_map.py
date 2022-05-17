"""
author: Ditlev Frickmann
email: reimer.frickmann@gmail.com

This script will create a segmentation map for all fits files beneath the specified root directory.

"""

from sys import argv
from shutil import move
from os.path import isdir

from MTLib import files
from MTLib.fitstools import create_segmentation_map

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

        '''
        We're going to create a mask of all the sources in the fits file.
        We use SExtractor to create a segmentation map, and then we create a source mask from that.
        '''
        segmentation_map = create_segmentation_map(fits_file)
        segmentation_map_out = segmentation_map.replace('Data','Output')
        print(f'Moving segmentation map to {segmentation_map_out}')
        move(segmentation_map, segmentation_map_out)


if __name__ == '__main__':
    main()