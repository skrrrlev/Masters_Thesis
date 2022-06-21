"""
author: Ditlev Frickmann
email: reimer.frickmann@gmail.com

This script will create a sigma image for the sourcefile <item> of a .ini file output tab.
Takes either a .ini file or a directory as argument.
    - If a .ini file is given, the script will only be applied to it.
    - If a directory is given, the script will be applied to all .ini files in and below that directory.
"""

from sys import argv
import numpy as np
from astropy.io import fits
from scipy.ndimage import binary_dilation

from MTLib.files import MapPipelineManager as MPP
from MTLib import PATH
from MTLib.fitstools import create_sigma_image_from_noise_distribution, create_from

SOURCE_MASK_PIXEL_DILATION = 7 # pixels


def main():

    ini_files = MPP.read_input(argv)
    ini_files = MPP.get_output_files(ini_files)

    for ini_file in ini_files:
        with MPP(ini_file) as ini:
            tab = 'DEFAULT'
            name = ini.get_item_from(tab, item='name')
            
            '''Make sure that the dead mask file exists and load in the mask'''
            if not ini.contains(tab,item='dead_mask_file'):
                print(f'Skipped creating sigma image for {name}. Could not find dead pixel mask.')
                continue
            dead_mask_file = ini.get_item_from(tab,item='dead_mask_file')
            with fits.open(dead_mask_file) as dead_hdul:
                dead_mask: np.ndarray = np.array(dead_hdul[1].data == 1, dtype=bool)

            '''
            Next, we're making sure that the segmentation map exists. 
            If it does, we create an initial source mask from the map.
            To avoid the few pixels included around the edge of sources, that are affected by the sources, we dilute the mask a little bit.
            '''
            if not ini.contains(tab,item='segmentation_map_file'):
                print(f'Skipped creating sigma image for {tab}. Could not find segmentation map.')
                continue
            segmentation_map_file = ini.get_item_from(tab,item='segmentation_map_file')
            with fits.open(segmentation_map_file) as sm_hdul:
                source_mask: np.ndarray = np.array(sm_hdul[1].data != 0, dtype=bool)
            source_mask = binary_dilation(source_mask, iterations=SOURCE_MASK_PIXEL_DILATION)

            '''
            Next, we create the mask of noise pixels, by considering all the pixels
            that are neither included in the dead pixel mask or the source mask.
            Then we save the mask, so that we can check if everything went right.
            '''
            noise_pixel_mask = np.logical_not(np.logical_or(dead_mask, source_mask))
            source_file = ini.get_item_from(tab,item='source_file')
            noise_mask_file = f"{ini.get_item_from(tab,item='path')}/{name}_noise_mask.fits"
            ini.add_item_to(tab, item='noise_mask_file', value=noise_mask_file)
            with create_from(from_file=source_file, to_file=noise_mask_file,index=0) as hdu:
                hdu.header['EXTVER'] = 'NOISEMASK'
                data = np.zeros(hdu.data.shape)
                data[noise_pixel_mask] = 1
                hdu.data = data

            '''Now that we have the pixel noise mask, we can create the sigma map.'''       
            sigma_file = create_sigma_image_from_noise_distribution(
                fits_file=source_file,
                noise_pixel_mask_file=noise_mask_file,
                dead_pixel_mask_file=dead_mask_file,
                noise_distribution_plot_file=PATH.FIGURES.value + 'noise_distribution/' + name + '_noise_distribution.png'
            )
            ini.add_item_to(tab, item='sigma_file', value=sigma_file)


if __name__ == '__main__':
    main()