"""
author: Ditlev Frickmann
email: reimer.frickmann@gmail.com

This script will create a sigma image for all fits files beneath the specified root directory if they match an input pattern.

"""

from sys import argv
from os.path import isdir, isfile
from os import listdir
import configparser

from MTLib import files
from MTLib.files import MapPipelineManager as MPP
from MTLib.fitstools import create_mask_from_segmentation_map


def setup():
    root = argv[1]
    if not isdir(root):
        raise ValueError('Root must be a directory.')
    return root

def main():

    ini_files = MPP.read_input(argv)
    ini_files = MPP.get_output_files(ini_files)

    for ini_file in ini_files:
        with MPP(ini_file) as ini:
            
            # Check whether some segments should be excluded from the map
            # as defined by the 'exclude_from_mask' keyword in the .ini file.
            exclude_from_mask = ini.get_item_from('DEFAULT',item='exclude_from_mask')
            exclude_list = [int(val) for val in exclude_from_mask.split(',') if val]
            
            # If the .ini contains the right ascension and declination, use those coordinates, to create the mask.
            # Otherwise, use the centre of the cutout.
            if ini.contains('DEFAULT','ra') and ini.contains('DEFAULT','dec'):
                right_ascension = float(ini.get_item_from('DEFAULT','ra'))
                declination = float(ini.get_item_from('DEFAULT','dec'))
            else:
                right_ascension = None
                declination = None

            segmentation_map_file = ini.get_item_from('galaxy',item='galaxy_segmentation_map_file')
            
            mask_file = create_mask_from_segmentation_map(
                segmentation_map=segmentation_map_file,
                ra=right_ascension,
                dec=declination,
                exclude_list=exclude_list
            )
            ini.add_item_to('galaxy',item='mask_file', value=mask_file)

if __name__ == '__main__':
    main()