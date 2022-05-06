"""
author: Ditlev Frickmann
email: reimer.frickmann@gmail.com

This script will apply crops to all fits files beneath the specified root directory.

"""

from sys import argv
from os.path import isdir

from MTLib import files
from MTLib.fitstools import crop_using_region

def setup():
    root = argv[1]
    if not isdir(root):
        raise ValueError('Root must be a directory.')
    return root

def main():
    root = setup()
    
    fits_files = files.walker(root,extension='.fits')
    reg_files = files.walker(root,extension='.reg')

    for fits_file in fits_files:
        file_name = files.extract_filename(fits_file)
        fits_file_path = files.extract_path(fits_file)
        output_dir = fits_file_path + 'cutout_' + file_name
        for reg_file in reg_files:
            if not file_name in reg_file:
                continue
            crop_using_region(
                fits_path=fits_file,
                region_path=reg_file,
                output_folder=output_dir,
                output_name=file_name
            )

if __name__ == '__main__':
    main()