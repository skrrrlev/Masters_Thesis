"""

"""
from sys import argv
from astropy.io import fits
import numpy as np

from MTLib.galfit import create_input_file
from MTLib.files import MapPipelineManager as MPP
from MTLib.fitstools import source_characteristics as SC

def main():
    print("\nCreating galfit input files...")
    ini_files = MPP.read_input(argv)
    ini_files = MPP.get_output_files(ini_files)

    for ini_file in ini_files:
        with MPP(ini_file) as ini:
            try:
                tab = 'galaxy'
                
                """ Set path and filename for the input file """
                name = ini.get_item_from(tab,item='name')
                path = ini.get_item_from(tab,item='path')
                input_file = f'{path}/{name}.input'

                """ Set the parameters for the input file """
                file_to_fit = ini.get_item_from(tab,item='galaxy_source_file')
                with fits.open(file_to_fit) as hdul:
                    hdu: fits.hdu.image._ImageBaseHDU = hdul[0]
                    N,M = hdu.data.shape # number of rows, number of columns
                    x_scale = float(hdu.header['CD2_2'])*3600
                    y_scale = float(hdu.header['CD1_1'])*3600
                    zero_point = hdu.header['ABMAGZP']
                    # zero-point: https://www.stsci.edu/hst/instrumentation/acs/data-analysis/zeropoints

                parameters = {
                    'A': [f'{file_to_fit}[0]'], # A: file to fit
                    'B': [f'{path}/{name}_output.fits'], # B: output file
                    'C': [f"{ini.get_item_from(tab,item='galaxy_sigma_file')}[1]"], # C: weight file
                    'D': [f"{ini.get_item_from(tab,item='psf_file')}[0]"], # D: psf file
                    'F': [f"{ini.get_item_from(tab,item='mask_file')}[0]"], # F: bad pixel mask
                    'H': ['1', f'{M}', '1', f'{N}'], # image region: xmin, xmax, ymin, ymax
                    'J': [f'{zero_point:.4f}'], # magnitude zero-point
                    'K': [f'{x_scale:.4f}', f'{y_scale:.4f}'] # plate scale in arcseconds per pixel

                }

                """ Add the items to fit to the image """
                # Get the approximate brightness of the source:
                brightness_mask_file=f'{path}/galaxy/brightness_mask.fits'
                approx_brightness = SC.estimate_brightness(
                    fits_file=file_to_fit,
                    sigma_file=ini.get_item_from(tab,item='galaxy_sigma_file'),
                    mask_file=ini.get_item_from(tab,item='mask_file'),
                    zero_point=zero_point,
                    brightness_mask_file=brightness_mask_file
                )

                # Get the position of the source in the image:
                source_col, source_row = SC.get_pixel_in_map_from_brightest_pixel(
                    map=file_to_fit,
                    extension=0,
                    brightness_mask=brightness_mask_file
                )

                # Get the approximate effective radius of the source in pixels
                approx_effective_radius = SC.estimate_half_light_radius(
                    map=file_to_fit,
                    extension=0,
                    brightness_mask=brightness_mask_file
                )

                sersic = {
                    '0': ['sersic'], # object type
                    '1': [f'{source_col}',f'{source_row}','1','1', '# x, y [pix]'], # position of source
                    '3': [f'{approx_brightness:.3f}','1','# mag'], # integrated magnitude
                    '4': [f'{approx_effective_radius:.3f}','1', '# re'], # effective radius
                    '5': ['2','1', '# n'], # sersic index
                    '9': ['0.5','1', '# AR'], # axis ratio
                    '10': ['0','1','# PA'], # position angle
                    'Z': ['0'] # output option (0=residual)
                }

                sky_bg = {
                    '0': ['sky'],
                    '1': ['0','1'],
                    '2': ['0.0000','0'],
                    '3': ['0.0000','0'],
                    'Z': ['0']
                } 

                """ Create the galfit input file"""
                create_input_file(input_file, parameters, [sersic,sky_bg])
                ini.add_item_to('DEFAULT',item='galfit_input',value=input_file)
            except Exception as e:
                print(f"Failed to create input file: {e}")




if __name__ == '__main__':
    main()