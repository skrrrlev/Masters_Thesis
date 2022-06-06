from shutil import copyfile
from os.path import isfile
from os import listdir, remove, getcwd, chdir

from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord, FK5, ICRS

from subprocess import Popen, PIPE
from numpy import isnan, round

from ..files import extract_filename, extract_path
from .. import PATH

class SegmentationMapSourceValueError(Exception):
    '''Raised when the value in the segmentation map where the source is supposed to be is zero.'''

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

        segmentation_file_name = fits_file_name+'_segmentation'
        new_content = content.replace('<segmentation>',segmentation_file_name)
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
    segmentation_file = fits_file_dir + segmentation_file_name + '.fits'
    if isfile(segmentation_file):
        remove(segmentation_file)
    print(f'Saving {segmentation_file}')
    chdir(fits_file_dir)
    try:
        process = Popen(['sex', fits_file_name+'.fits'], stdout=PIPE, stderr=PIPE)
        process.communicate()
    except Exception as e:
        chdir(out_of_scope_directory)
        raise e(f'Failed to save {segmentation_file_name}.')

    '''Remove unessesary files and return the path to the segmentation map'''
    chdir(out_of_scope_directory)
    new_items_in_dir = [item for item in listdir(fits_file_dir) if item not in items_in_dir and 'segmentation' not in item]
    for item in new_items_in_dir:
        remove(fits_file_dir+item)
    return segmentation_file

def create_mask_from_segmentation_map(segmentation_map: str, ra: float = None, dec: float = None, frame=FK5, exclude_list:"list[int]"=[]):
    '''
    Using a segmentation map and the coordinates of the source, create a mask for galfit.
    If coordinates are not specified, will use the centre of the map.
    '''

    '''Load the HDU and the WCS'''
    hdu = fits.open(segmentation_map)[0]
    wcs = WCS(hdu.header)

    '''Get the pixels at the centre of the source'''
    N,M = hdu.data.shape
    if ra is None or dec is None:
        x = int(N/2)
        y = int(M/2)
    else:
        sky = SkyCoord(ra,dec,frame=frame, unit='deg')
        xp, yp = wcs.world_to_pixel(sky)
        if isnan(xp) or isnan(yp):
            raise ValueError('Input coordinates is not within image.')
        x = int(round(xp,0))
        y = int(round(M-yp,0))

    '''Get the value of the segmentation map at the source'''
    value_at_source = hdu.data[x,y]
    if not value_at_source:
        raise SegmentationMapSourceValueError(f'Source value at ({x},{y}) in segmentation map was zero.')

    '''Remove the source from the mask'''
    hdu.data[hdu.data == value_at_source] = 0

    '''Remove the sources in the exclude list'''
    for value in exclude_list:
        hdu.data[hdu.data == value] = 0

    '''Save the mask'''
    mask_file_name = segmentation_map.replace('segmentation','mask')
    print(f'Saving {mask_file_name}')
    hdu.writeto(mask_file_name, overwrite=True)

if __name__ == '__main__':
    a = create_segmentation_map(PATH.OUTPUTMAPS.value + 'HST/source13/source13_F160w/source13_F160w_galaxy.fits')
    #print(a)
    create_mask_from_segmentation_map('Output/Maps/HST/source13/source13_F160w/source13_F160w_galaxy_segmentation.fits',150.07597,2.2118269,frame=ICRS)
    #create_mask_from_segmentation_map('Output/Maps/HST/source13/source13_F814w/source13_F814w_galaxy_segmentation.fits',2,2.2118269)