from re import findall
from os import makedirs
from os.path import exists, isdir
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
import configparser
from .create_from import create_from
from ..files import MapPipelineManager as MPP

def crop_using_region(fits_path: str, region_path: str, output_folder: str, output_name: str, extension:int=0) -> "list[str]":
    '''Use a ds9 region file in image units to cut out a crop of a fits file, and save it as a new fits file with a correct header.'''

    with open(region_path,'r') as f:
        content = f.read().lower().split('\n')

    crops: "dict[str,dict[str,float]]" = {}
    for line in content:
        matches = findall(r'(\d+(?:\.\d+)?(?:e[+-]?\d+)?)\S', line)
        if len(matches) >= 5:
            try:
                region_name = findall(r'text=\{(\S+)\}', line)[0]
            except IndexError: # A name was not defined for the region
                continue
            crops[region_name] = {
                'position' : (float(matches[0])-1,float(matches[1])-1), # -1 to account for zero-based index.
                'size' : (float(matches[2]),float(matches[3])),
                'angle' : float(matches[4]),
            }
    
        for invalid_type in ['icrs','fk5','fk4']:
            if invalid_type in line:
                raise ValueError('Region map is not saved in image units.')
    
    files = []
    for region_name in crops:
        crop = crops[region_name]
        hdu = fits.open(fits_path)[extension]
        wcs = WCS(hdu.header)

        # Make the cutout, including the WCS
        cutout = Cutout2D(hdu.data, position=crop['position'], size=crop['size'], wcs=wcs)

        # Put the cutout image in the FITS HDU
        hdu.data = cutout.data

        # Update the FITS header with the cutout WCS
        hdu.header.update(cutout.wcs.to_header())

        # Write the cutout to a new FITS file
        crop_output_folder = output_folder + '/'*(not output_folder.endswith('/')) + region_name + '/'
        if not exists(crop_output_folder):
            makedirs(crop_output_folder)
        cutout_filename = crop_output_folder + output_name + '.fits'
        print(f'Saving {cutout_filename}')
        hdu.writeto(cutout_filename, overwrite=True)
        files.append(cutout_filename)
    return files

def get_regions(region_file:str) -> "list[str]":
    '''Get the names of the regions defined in a .reg region file.'''
    with open(region_file,'r') as f:
        content = f.read().lower().split('\n')
    
    regions: list[str] = []
    for line in content:
        matches = findall(r'(\d+(?:\.\d+)?(?:e[+-]?\d+)?)\S', line)
        if len(matches) >= 5:
            try:
                region_name = findall(r'text=\{(\S+)\}', line)[0]
                if region_name == 'out':
                    continue
                regions.append(region_name)
            except IndexError: # A name was not defined for the region
                continue

    return regions

def crop_using_ini_tab(ini_tab:configparser.SectionProxy, file_key: str, name:str, extension:int=0) -> None:
    output_path = ini_tab["path"]
    region_path = ini_tab["region"]

    with open(region_path,'r') as f:
        content = f.read().lower().split('\n')

    crops: "dict[str,dict[str,float]]" = {}
    for line in content:
        matches = findall(r'(\d+(?:\.\d+)?(?:e[+-]?\d+)?)\S', line)
        if len(matches) >= 5:
            try:
                region_name = findall(r'text=\{(\S+)\}', line)[0]
                if region_name == 'out':
                    continue
            except IndexError: # A name was not defined for the region
                continue
            crops[region_name] = {
                'position' : (float(matches[0])-1,float(matches[1])-1), # -1 to account for zero-based index.
                'size' : (float(matches[2]),float(matches[3])),
                'angle' : float(matches[4]),
            }
    
        for invalid_type in ['icrs','fk5','fk4']:
            if invalid_type in line:
                raise ValueError('Region map is not saved in image units.')

    files = []
    for region_name in crops:
        crop = crops[region_name]
        crop_output_folder = f'{output_path}/{region_name}'
        if not isdir(crop_output_folder):
            makedirs(crop_output_folder)
        output_file = f'{crop_output_folder}/{name}.fits'
        with create_from(ini_tab["region"][file_key],output_file, index=extension) as hdu:
            wcs = WCS(hdu.header)
            # Make the cutout, including the WCS
            cutout = Cutout2D(hdu.data, position=crop['position'], size=crop['size'], wcs=wcs)
            # Put the cutout image in the FITS HDU
            hdu.data = cutout.data
            # Update the FITS header with the cutout WCS
            hdu.header.update(cutout.wcs.to_header())


        files.append(output_file)
    return files

def crop_using_MPP(ini: MPP, file_key:str, output_name:str, extension:int=0):
    '''Apply crops defined by the region file in a .ini MPP, to the file given by the <file_key>.'''
    region_path = ini.get_item_from('DEFAULT',item='region_file')

    with open(region_path,'r') as f:
        content = f.read().lower().split('\n')

    crops: "dict[str,dict[str,float]]" = {}
    for line in content:
        matches = findall(r'(\d+(?:\.\d+)?(?:e[+-]?\d+)?)\S', line)
        if len(matches) >= 5:
            try:
                region_name = findall(r'text=\{(\S+)\}', line)[0]
                if region_name == 'out':
                    continue
            except IndexError: # A name was not defined for the region
                continue
            crops[region_name] = {
                'position' : (float(matches[0])-1,float(matches[1])-1), # -1 to account for zero-based index.
                'size' : (float(matches[2]),float(matches[3])),
                'angle' : float(matches[4]),
            }
    
        for invalid_type in ['icrs','fk5','fk4']:
            if invalid_type in line:
                raise ValueError('Region map is not saved in image units.')

    for region_name in crops:
        crop = crops[region_name]
        crop_output_folder = ini.get_item_from(region_name, item=f'{region_name}_path')
        if not isdir(crop_output_folder):
            makedirs(crop_output_folder)
        
        input_file = ini.get_item_from('DEFAULT',item=file_key)
        output_file = f'{crop_output_folder}/{output_name}.fits'
        with create_from(from_file=input_file, to_file=output_file, index=extension) as hdu:
            wcs = WCS(hdu.header)
            # Make the cutout, including the WCS
            cutout = Cutout2D(hdu.data, position=crop['position'], size=crop['size'], wcs=wcs)
            # Put the cutout image in the FITS HDU
            hdu.data = cutout.data
            # Update the FITS header with the cutout WCS
            hdu.header.update(cutout.wcs.to_header())
        ini.add_item_to(region_name, item=f'{region_name}_{file_key}', value=output_file)



if __name__ == '__main__':
    '''
    crop_using_region(
        fits_path='Data/Maps/HST/source13/source13_F160w.fits',
        region_path='Data/Maps/HST/source13/regions/source13_F160w.reg',
        output_folder='Data/Maps/HST/source13/cutout_source13_F160w/',
        output_name='source13_F160w'
    )
    '''

    with MPP('Output/Maps/s14_F140w/s14_F140w.ini') as ini:
        crop_using_MPP(ini,'source_file','test',extension=0)
