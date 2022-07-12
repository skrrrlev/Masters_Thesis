"""

"""
from sys import argv

from MTLib.files import MapPipelineManager as MPP
from MTLib.galfit import run_galfit


def main():
    print('\nRunning galfit...')
    ini_files = MPP.read_input(argv)
    ini_files = MPP.get_output_files(ini_files)

    failed = []
    for ini_file in ini_files:
        with MPP(ini_file) as ini:
            input_file = ini.get_item_from('DEFAULT','galfit_input')
            path = ini.get_item_from('DEFAULT','path')
            name = ini.get_item_from('DEFAULT','name')
            output_file = f'{path}/{name}_output.fits'
            log_file = f'{ini.get_item_from("DEFAULT",item="path")}/{ini.get_item_from("DEFAULT",item="name")}_log.txt'
            try:
                run_galfit(input_file, output_file=output_file, log_file=log_file)
            except Exception as e:
                print(e)
                failed.append(ini.get_item_from('DEFAULT',item='name'))

    print(f'Failed to fit adequately: {", ".join(failed)}')

if __name__ == '__main__':
    main()