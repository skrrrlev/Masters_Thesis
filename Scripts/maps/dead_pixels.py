"""
author: Ditlev Frickmann
email: reimer.frickmann@gmail.com

This script will create a dead pixel mask from the weight image in the second extension of the fits file.
Takes either a .ini file or a directory as argument.
    - If a .ini file is given, the script will only be applied to it.
    - If a directory is given, the script will be applied to all .ini files in and below that directory.

"""

from sys import argv
import numpy as np
from MTLib.files import MapPipelineManager as MPP
from MTLib.fitstools import create_from
from os import remove
from os.path import isfile

def main():
    ini_files = MPP.read_input(argv)
    ini_files = MPP.get_output_files(ini_files)

    for ini_file in ini_files:
        with MPP(ini_file) as ini:
            tab='DEFAULT'
            path = ini.get_item_from(tab,item="path")
            name = ini.get_item_from(tab,item="name")
            dead_mask = f'{path}/{name}_dead_mask.fits'
            if isfile(dead_mask):
                remove(dead_mask)
            ini.add_item_to(tab,item='dead_mask_file',value=dead_mask)

            with create_from(from_file=ini.get_item_from(tab,item="source_file"), to_file=dead_mask, index=1) as wht_hdu:
                dead_pixel_mask = wht_hdu.data == 0

                data = np.zeros(wht_hdu.data.shape)
                data[dead_pixel_mask] = 1

                wht_hdu.data = data
                wht_hdu.header['EXTVER'] = 'DEADMASK'


if __name__ == '__main__':
    main()