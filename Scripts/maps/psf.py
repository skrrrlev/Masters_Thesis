from sys import argv
from os.path import isdir
from os import listdir
from multiprocessing import Pool, cpu_count, freeze_support

from MTLib import files
from MTLib.files import MapPipelineManager as MPP
from MTLib.fitstools import create_psf

def setup():
    try:
        root = argv[1]
    except IndexError:
        raise ValueError('Please specify a root directory')
    if not isdir(root):
        raise ValueError('Root must be a directory.')
    return root
    

def main():
    print('Creating PSFs...')
    ini_files = MPP.read_input(argv)
    ini_files = MPP.get_output_files(ini_files)
    
    for ini_file in ini_files:
        with MPP(ini_file) as ini:
            name = ini.get_item_from('DEFAULT',item='name')
            if not 'psf' in ini.config:
                print(f'The region "psf" was not defined for {name}. Skipping creation of psf.')
                continue
            psf_file = f"{ini.get_item_from('galaxy',item='galaxy_path')}/{name}_psf.fits"
            create_psf(
                fits_file=ini.get_item_from('psf','psf_source_file'),
                sigma_file=ini.get_item_from('psf','psf_sigma_file'),
                region_file=ini.get_item_from('DEFAULT','region_file'),
                smooth_size=3,
                output_file=psf_file
            )
            ini.add_item_to('galaxy',item='psf_file',value=psf_file)
        

if __name__ == '__main__':
    main()