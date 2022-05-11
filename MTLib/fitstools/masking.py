
from shutil import copyfile
from os.path import isfile
from os import listdir, remove, getcwd, chdir

from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord, FK5, ICRS

from subprocess import Popen, PIPE
from numpy import isnan, round

from MTLib.files import extract_filename, extract_path
from MTLib import PATH

def create_segmentation_map(fits_file: str) -> str:
    '''Use SExtractor to create a segmentation map for a fits file.'''

    out_of_scope_directory = getcwd()

    '''Get the directory and name of the fits file'''
    fits_file_dir = extract_path(fits_file)
    fits_file_name = extract_filename(fits_file)

    items_in_dir = listdir(fits_file_dir)

    '''Get the SExtractor configuration files'''
    configuration_dir = extract_path(__file__) + 'masktools/'
    sex_file = configuration_dir + 'default.sex'
    param_file = configuration_dir + 'default.param'

    try:

        '''Move and adapt the .sex file'''
        with open(sex_file, 'r') as sf:
            content = sf.read()
        
        new_sex_file = fits_file_dir + 'default.sex'
        if isfile(new_sex_file):
            remove(new_sex_file)

        new_content = content.replace('<segmentation>',fits_file_name+'_segmentation')
        with open(new_sex_file,'w') as nsf:
            nsf.write(new_content)
        
        '''Copy the params file'''
        copyfile(param_file, fits_file_dir+'default.param')

    except Exception as e:

        '''If something goes wrong, delete the new files and return to the orginal working directory before raising the Exception.'''
        new_items_in_dir = [item for item in listdir(fits_file_dir) if item not in items_in_dir]
        for item in new_items_in_dir:
            remove(fits_file_dir+item)
        
        raise e

    '''Change directory and run SExtractor'''
    chdir(fits_file_dir)
    print(getcwd())
    try:
        process = Popen(['sex', fits_file_name+'.fits'], stdout=PIPE, stderr=PIPE)
        process.communicate()
    except Exception as e:
        chdir(out_of_scope_directory)
        raise e

    '''Remove unessesary files and return the path to the segmentation map'''
    new_items_in_dir = [item for item in listdir(fits_file_dir) if item not in items_in_dir and 'segmentation' not in item]
    for item in new_items_in_dir:
        remove(item)
    chdir(out_of_scope_directory)
    segmentation_file = fits_file_dir + fits_file_name+'_segmentation.fits'
    print(f'Saving {segmentation_file}')
    return segmentation_file

def create_mask_from_segmentation_map(segmentation_map: str, ra: float, dec: float, frame=FK5):

    '''Load the HDU and the WCS'''
    hdu = fits.open(segmentation_map)[0]
    wcs = WCS(hdu.header)

    '''Get the pixels at the centre of the source'''
    sky = SkyCoord(ra,dec,frame=frame, unit='deg')
    x, y = wcs.world_to_pixel(sky)
    if isnan(x) or isnan(y):
        raise ValueError('Input coordinates is not within image.')

    '''Get the value of the segmentation map at the source'''
    N,M = hdu.data.shape
    value_at_source = hdu.data[int(round(x,0)),int(round(M-y,0))]
    
    '''Remove the source from the mask'''
    hdu.data[hdu.data == value_at_source] = 0

    '''Save the mask'''
    mask_file_name = extract_path(segmentation_map) + extract_filename(segmentation_map).replace('segmentation','mask') + '.fits'
    print(f'Saving {mask_file_name}')
    hdu.writeto(mask_file_name, overwrite=True)



if __name__ == '__main__':
    a = create_segmentation_map(PATH.OUTPUTMAPS.value + 'HST/source13/source13_F160w/source13_F160w_galaxy.fits')
    #print(a)
    create_mask_from_segmentation_map('Output/Maps/HST/source13/source13_F160w/source13_F160w_galaxy_segmentation.fits',150.07597,2.2118269,frame=ICRS)
    #create_mask_from_segmentation_map('Output/Maps/HST/source13/source13_F814w/source13_F814w_galaxy_segmentation.fits',2,2.2118269)