"""
author: Ditlev Frickmann
email: reimer.frickmann@gmail.com

This script will apply crops to all fits files beneath the specified root directory if they have a corresponding region file.

"""

from sys import argv
from os.path import isdir, isfile
from os import makedirs
from shutil import rmtree

from matplotlib.pyplot import get

from MTLib.files import MapPipelineManager as MPP
from MTLib import files
from MTLib.fitstools import crop_using_region, get_regions
from MTLib.fitstools.croptools import crop_using_MPP, crop_using_ini_tab


def main():
    ini_files = MPP.read_input(argv)
    ini_files = MPP.get_output_files(ini_files)

    for ini_file in ini_files:
        with MPP(ini_file) as ini:
            tab = 'DEFAULT'
            path = ini.get_item_from(tab, item="path")
            name = ini.get_item_from(tab, item="name")
            
            if not ini.contains(tab, item="region_file"):
                print(f'Skipped creating crops for {name}. Region file did not exist.')
                continue
            region_file = ini.get_item_from(tab, item="region_file")
            regions = get_regions(region_file)
            
            for region in regions:
                ini.add_tab(region)

                region_path = f'{path}/{region}'
                if isdir(region_path):
                    rmtree(region_path)
                makedirs(region_path)
                ini.add_item_to(region, item=f'{region}_path', value=region_path)

            crop_using_MPP(ini, file_key='source_file', output_name=f'{name}', extension=0)
            crop_using_MPP(ini, file_key='dead_mask_file', output_name=f'{name}_dead_mask', extension=1)
            crop_using_MPP(ini, file_key='segmentation_map_file', output_name=f'{name}_segmentation_map', extension=1)
            crop_using_MPP(ini, file_key='sigma_file', output_name=f'{name}_sigma', extension=1)
                

if __name__ == '__main__':
    main()