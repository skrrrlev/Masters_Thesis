from sys import argv
from os.path import isdir
from os import listdir
from multiprocessing import Pool, cpu_count, freeze_support

from MTLib import files
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
    root = setup()

    # get fits files in root and below root
    fits_files = [file for file in files.walker(root,extension='.fits')]

    psf_files, sigma_files, region_files = [],[],[]
    for fits_file in fits_files:
        # get name of current source from the fits file
        file_name = files.extract_filename(fits_file)
        fits_dir = files.extract_path(fits_file)
        output_dir = fits_dir.replace('Data','Output')

        temp_psf_files = [output_dir+item+'/'+file_name+'.fits' for item in listdir(output_dir) if isdir(output_dir+item+'/') and 'psf' in item]
        temp_sigma_files = [item.replace('.fits','_sigma.fits') for item in temp_psf_files]
        temp_region_files = [fits_dir+f'regions/{file_name}.reg' for _ in range(len(temp_psf_files))]

        psf_files += temp_psf_files
        sigma_files += temp_sigma_files
        region_files += temp_region_files
    
    for args in zip(psf_files,sigma_files,region_files):
        try:
            create_psf(*args)
        except Exception as e:
            print(f'Error generating PSF for "{args[0]}": {e}')
        

if __name__ == '__main__':
    main()