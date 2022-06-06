"""
author: Ditlev Frickmann
email: reimer.frickmann@gmail.com

This script will create a sigma image for all fits files beneath the specified root directory if they match an input pattern.

"""

from sys import argv
from os.path import isdir, isfile
from os import listdir
import configparser

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

        exclude_list = []
        ini_file = output_dir+file_name+'.ini'
        if isfile(ini_file):
            config = configparser.ConfigParser()
            config.read(ini_file)
            if 'MASK' in config:
                if 'exclude' in config['MASK'].keys():
                    exclude = config['MASK']['exclude']
                    exclude_list = [int(val) for val in exclude.split(',')]
            else:
                config['MASK'] = {}
            with open(ini_file, 'w') as config_file:
                config.write(config_file)
        print(f'Exclude in ini: {exclude_list}')

        cutout_directories = [output_dir+item+'/' for item in listdir(output_dir) if isdir(output_dir+item+'/')]

        for subdir in cutout_directories:
        
            try:
                segmentation_map = subdir + file_name + '_segmentation.fits'
                create_mask_from_segmentation_map(segmentation_map, exclude_list)
            except Exception as e:
                print(f'Failed to create mask for {segmentation_map}: {e}')

if __name__ == '__main__':
    main()