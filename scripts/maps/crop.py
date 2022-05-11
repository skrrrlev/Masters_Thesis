"""
author: Ditlev Frickmann
email: reimer.frickmann@gmail.com

This script will apply crops to all fits files beneath the specified root directory if they have a corresponding region file.

"""

from sys import argv
from os.path import isdir

from MTLib import files
from MTLib.fitstools import crop_using_region
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
    # get region files in root and below
    reg_files = files.walker(root,extension='.reg')

    for fits_file in fits_files:
        # get name of current source from the fits file
        file_name = files.extract_filename(fits_file)
        # define the output folder
        fits_file_path = files.extract_path(fits_file)
        output_dir = fits_file_path.replace('Data','Output') + file_name
        # find the corresponding region file
        for reg_file in reg_files:
            if not file_name in reg_file:
                continue
            # cut out the regions described in the region file.
            crop_using_region(
                fits_path=fits_file,
                region_path=reg_file,
                output_folder=output_dir,
                output_name=file_name
            )

if __name__ == '__main__':
    main()