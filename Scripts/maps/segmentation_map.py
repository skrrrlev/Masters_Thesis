"""
author: Ditlev Frickmann
email: reimer.frickmann@gmail.com

This script will create a segmentation map for all fits files beneath the specified root directory.
Takes either a .ini file or a directory as argument.
    - If a .ini file is given, the script will only be applied to it.
    - If a directory is given, the script will be applied to all .ini files in and below that directory.
"""

from sys import argv
from MTLib.files import MapPipelineManager as MPP
from MTLib.fitstools import create_segmentation_map


def main():
    print("\nCreating segmentation maps...")
    ini_files = MPP.read_input(argv)
    ini_files = MPP.get_output_files(ini_files)

    failures = []
    for ini_file in ini_files:
        with MPP(ini_file) as ini:
            tab='DEFAULT'
            name = ini.get_item_from(tab,'name')
            source_file = ini.get_item_from(tab, 'source_file')
            try:
                segmentation_map = create_segmentation_map(source_file)
                ini.add_item_to(tab,item='segmentation_map_file',value=segmentation_map)
            except Exception as e:
                failures.append(f"{name}: {e}")
    
    if failures:
        print("Unable to create segmentation map for the following:")
        for e in failures:
            print(e)


if __name__ == '__main__':
    main()