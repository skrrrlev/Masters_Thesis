'''
author: Ditlev Frickmann

Use the Pipeline manager for the data processing of the maps to initialize the output structure.
Takes either a .ini file or a directory as argument.
    - If a .ini file is given, the initialization will only be applied to it.
    - If a directory is given, the initialization will be applied to all .ini files in and below that directory.
'''
from sys import argv
from MTLib.files import MapPipelineManager as MPP

def main():
    print("\nInitializing pipeline...")
    ini_files = MPP.read_input(argv)
    for ini_file in ini_files:
        MPP.initialize_from_file(ini_file)

if __name__ == '__main__':
    main()