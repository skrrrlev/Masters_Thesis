from re import findall
from os import makedirs
from os.path import exists
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import Cutout2D


def crop_using_region(fits_path: str, region_path: str, output_folder: str, output_name: str):
    '''Use a ds9 region file in image units to cut out a crop of a fits file, and save it as a new fits file with a correct header.'''
    if not exists(output_folder):
        makedirs(output_folder)

    with open(region_path,'r') as f:
        content = f.read().lower().split('\n')

    crops: "dict[str,dict[str,float]]" = {}
    for line in content:
        matches = findall(r'(\d+(?:\.\d+)?(?:e[+-]?\d+)?)\S', line)
        if len(matches) == 5:
            region_name = findall(r'text=\{(\S+)\}', line)[0]
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
        hdu = fits.open(fits_path)[0]
        wcs = WCS(hdu.header)

        # Make the cutout, including the WCS
        cutout = Cutout2D(hdu.data, position=crop['position'], size=crop['size'], wcs=wcs)

        # Put the cutout image in the FITS HDU
        hdu.data = cutout.data

        # Update the FITS header with the cutout WCS
        hdu.header.update(cutout.wcs.to_header())

        # Write the cutout to a new FITS file
        cutout_filename = output_folder + '/'*(not output_folder.endswith('/')) + output_name + '_' + region_name + '.fits'
        print(f'Saving {cutout_filename}')
        hdu.writeto(cutout_filename, overwrite=True)

if __name__ == '__main__':
    crop_using_region(
        fits_path='Data/Maps/HST/source13/source13_F160w.fits',
        region_path='Data/Maps/HST/source13/regions/source13_F160w.reg',
        output_folder='Data/Maps/HST/source13/cutout_source13_F160w/',
        output_name='source13_F160w'
    )