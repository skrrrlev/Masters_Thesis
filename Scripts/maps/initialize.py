
from sys import argv
import numpy as np
from os.path import isfile, isdir
from os import listdir, makedirs
from shutil import copyfile
from astropy.io import fits

from MTLib import files

class InitializationError(Exception):
    '''Error raised when initialisation fails.'''

def setup():
    root = argv[1]
    if not isdir(root):
        raise ValueError('Root must be a directory.')
    return root

def main():
    root = setup()
    
    print(f'Initializing pipeline for all fits file below {root}')
    # get fits files in root and below root
    fits_files = files.walker(root,extension='.fits')

    for fits_file in fits_files:
        '''get name and path of the current fits file'''
        file_name = files.extract_filename(fits_file)
        fits_file_path = files.extract_path(fits_file)

        
        output_path = fits_file_path.replace('Data','Output')
        if not isdir(output_path):
            makedirs(output_path)

        ini_list = [item for item in listdir(fits_file_path) if item.endswith('.ini')]

        if not len(ini_list) == 1:
            raise InitializationError(f'Expected 1 .ini file in {fits_file_path}, but found {len(ini_list)}')

        ini_file = ini_list[0]

        copyfile(fits_file_path+ini_file,output_path+file_name+'.ini')
        print(f'Initialized {fits_file}')

if __name__ == '__main__':
    main()